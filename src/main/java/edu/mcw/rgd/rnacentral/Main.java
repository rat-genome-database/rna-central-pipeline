package edu.mcw.rgd.rnacentral;

import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * @author mtutaj
 * @since 3/01/19
 */
public class Main {

    private DAO dao = new DAO();
    private String version;
    private String pipelineName;
    private int xdbKeyForRNACentral;
    private String refSeqMappingFile;
    private String rgdMappingFile;

    Logger log = LogManager.getLogger("status");
    Logger logMultimatch = LogManager.getLogger("multimatch");

    public static void main(String[] args) throws Exception {

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Main manager = (Main) (bf.getBean("manager"));

        try {
            manager.run();
        }catch (Exception e) {
            Utils.printStackTrace(e, manager.log);
            throw e;
        }
    }

    public void run() throws Exception {

        long startTime = System.currentTimeMillis();

        String msg = getVersion();
        log.info(msg);

        msg = dao.getConnectionInfo();
        log.info("   "+msg);

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        log.info("   started at "+sdt.format(new Date(startTime)));

        // download RefSeq file
        FileDownloader fd = new FileDownloader();
        fd.setExternalFile(getRefSeqMappingFile());
        fd.setLocalFile("data/refseq_mapping");
        fd.setUseCompression(true);
        fd.setPrependDateStamp(true);
        String localFile = fd.downloadNew();

        List<Integer> speciesTypeKeys = new ArrayList<>(SpeciesType.getSpeciesTypeKeys());
        speciesTypeKeys.removeIf(speciesTypeKey -> !SpeciesType.isSearchable(speciesTypeKey));

        // process every species in parallel
        CounterPool counters = new CounterPool();
        speciesTypeKeys.parallelStream().forEach( speciesTypeKey -> {
            try {
                int idCount = run(speciesTypeKey, localFile);
                counters.add(SpeciesType.getCommonName(speciesTypeKey), idCount);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        });

        // TODO: download ensembl file and process it as well

        log.info("");
        log.info("=== RNACentral id count");
        for( int speciesTypeKey: speciesTypeKeys ) {
            String speciesName = SpeciesType.getCommonName(speciesTypeKey);
            int count = counters.get(speciesName);
            if( count>0 ) {
                log.info(String.format("%12s - %7d", speciesName, count));
            }
        }

        log.info("");
        log.info("===    time elapsed: "+ Utils.formatElapsedTime(startTime, System.currentTimeMillis()));
        log.info("");
    }

    public int run(int speciesTypeKey, String localFile) throws Exception {

        // get taxon id for given species
        final String taxonId = Integer.toString(SpeciesType.getTaxonomicId(speciesTypeKey));
        CounterPool counters = new CounterPool();

        String species = SpeciesType.getCommonName(speciesTypeKey);
        log.debug("START for "+species);


        List<XdbId> idsIncoming = new ArrayList<>();

        // parse incoming files
        parseRefSeqFile(idsIncoming, localFile, taxonId, species, counters);
        if( speciesTypeKey==SpeciesType.RAT ) {
            parseRgdFile(idsIncoming, taxonId, species, counters);
        }


        // QC
        List<XdbId> idsInRgd = dao.getRNACentralIds(speciesTypeKey, getPipelineName(), getXdbKeyForRNACentral());
        log.debug("QC: get incoming "+getPipelineName()+" Ids for "+species);

        // determine to-be-inserted RNACentral ids
        log.debug("QC: determine to-be-inserted "+getPipelineName()+" Ids for "+species);
        Collection<XdbId> idsToBeInserted = CollectionUtils.subtract(idsIncoming, idsInRgd);

        // determine matching RNACentral ids
        log.debug("QC: determine matching "+getPipelineName()+" Ids for "+species);
        Collection<XdbId> idsMatching = CollectionUtils.intersection(idsIncoming, idsInRgd);
        idsMatching.retainAll(idsInRgd);

        // determine to-be-deleted RNACentral ids
        log.debug("QC: determine to-be-deleted "+getPipelineName()+" Ids for "+species);
        idsInRgd.removeAll(idsIncoming);
        Collection<XdbId> idsToBeDeleted = CollectionUtils.subtract(idsInRgd, idsIncoming);


        // loading
        if( !idsToBeInserted.isEmpty() ) {
            dao.insertXdbs(idsToBeInserted);
        }

        if( !idsToBeDeleted.isEmpty() ) {
            dao.deleteXdbIds(idsToBeDeleted);
        }

        if( !idsMatching.isEmpty() ) {
            dao.updateModificationDate(idsMatching);
        }


        synchronized(this) {
            log.info("===");
            log.info("summary for "+species);
            log.info("   RefSeq lines processed = " + Utils.formatThousands(counters.get("linesProcessedForSpecies")));

            int matchByRefSeq = counters.get("matchByRefSeq");
            if( matchByRefSeq!=0 ) {
                log.info("   match by acc id      = " + Utils.formatThousands(matchByRefSeq));
            }

            int noMatchByRefSeq = counters.get("noMatchByRefSeq");
            if( noMatchByRefSeq!=0 ) {
                log.info("   no match by acc id   = " + Utils.formatThousands(noMatchByRefSeq));
            }

            int multimatchByRefSeq = counters.get("noMatchByRefSeq");
            if( multimatchByRefSeq!=0 ) {
                log.info("   multimatch by acc id = " + Utils.formatThousands(multimatchByRefSeq));
            }


            int matchByRgdId = counters.get("matchByRgdId");
            if( matchByRgdId!=0 ) {
                log.info("   match by rgd id      = " + Utils.formatThousands(matchByRgdId));
            }

            int noMatchByRgdId = counters.get("noMatchByRgdId");
            if( noMatchByRgdId!=0 ) {
                log.info("   no match by rgd id   = " + Utils.formatThousands(noMatchByRgdId));
            }

            int matchByRgdIdDuplicates = counters.get("matchByRgdIdDuplicates");
            if( matchByRgdIdDuplicates!=0 ) {
                log.info("   match by rgd id dups = " + Utils.formatThousands(matchByRgdIdDuplicates));
            }


            if( idsToBeInserted.size() + idsToBeDeleted.size() + idsMatching.size() > 0 ) {
                log.info("");

                if (!idsToBeInserted.isEmpty()) {
                    log.info("  inserted " + getPipelineName() + " ids : " + Utils.formatThousands(idsToBeInserted.size()));
                }
                if (!idsToBeDeleted.isEmpty()) {
                    log.info("  deleted " + getPipelineName() + " ids : " + Utils.formatThousands(idsToBeDeleted.size()));
                }
                if (!idsMatching.isEmpty()) {
                    log.info("  matching " + getPipelineName() + " ids : " + Utils.formatThousands(idsMatching.size()));
                }
            }
        }

        return idsMatching.size() + idsToBeInserted.size() - idsToBeDeleted.size();
    }

    void parseRefSeqFile( List<XdbId> idsIncoming, String localFile, String taxonId, String species, CounterPool counters ) throws Exception {
        // file content
        //Tab-separated file with RNAcentral ids, corresponding external ids,
        //NCBI taxon ids, RNA types (according to INSDC classification),
        //        and gene names.
        // URS0000008E6C	REFSEQ	NR_113675	7	    rRNA
        // URS00008120E8	REFSEQ	NR_131204	9606	lncRNA	XACT

        BufferedReader in = Utils.openReader(localFile);
        String line;
        while( (line=in.readLine())!=null ) {
            String[] cols = line.split("[\\t]", -1);
            String rnaCentralId = cols[0];
            String tag = cols[1]; // must be 'REFSEQ'
            String accId = cols[2]; // RefSeq acc, must be 'NR_xxx..'
            String taxon = cols[3];
            String rnaType = cols[4];
            String geneSymbol = cols[5];

            if( !taxonId.equals(taxon) ) {
                continue;
            }
            if( !tag.equals("REFSEQ") ) {
                log.error("*** unexpected db tag: "+tag+";  was expecting:REFSEQ");
                continue;
            }
            counters.increment("linesProcessedForSpecies");

            List<Transcript> transcripts = dao.getTranscriptsByAccId(accId);
            if( transcripts.isEmpty() ) {

                List<Gene> genes = dao.getActiveGeneIdsForRefseqAcc(accId);
                if( genes.isEmpty() ) {
                    counters.increment("noMatchByRefSeq");
                    log.debug("-- no match for " + accId + " gene " + geneSymbol + "  species " + species);
                } else if( genes.size()==1 ) {
                    counters.increment("matchByRefSeq");

                    XdbId x = new XdbId();
                    x.setAccId(rnaCentralId);
                    x.setSrcPipeline(getPipelineName());
                    x.setRgdId(genes.get(0).getRgdId());
                    x.setXdbKey(getXdbKeyForRNACentral());
                    x.setCreationDate(new Date());
                    x.setModificationDate(new Date());
                    idsIncoming.add(x);
                } else {
                    counters.increment("multimatchByRefSeq");

                    String info = null;
                    for( Gene g: genes ) {
                        if( info==null ) {
                            info = g.getSymbol()+" (RGD:"+g.getRgdId()+")";
                        } else {
                            info += "  , "+g.getSymbol()+" (RGD:"+g.getRgdId()+")";
                        }
                    }
                    logMultimatch.debug(species+": "+accId+" matches multiple genes: "+info);
                }
            } else if( transcripts.size()==1 ) {
                counters.increment("matchByRefSeq");

                XdbId x = new XdbId();
                x.setAccId(rnaCentralId);
                x.setSrcPipeline(getPipelineName());
                x.setRgdId(transcripts.get(0).getGeneRgdId());
                x.setXdbKey(getXdbKeyForRNACentral());
                x.setCreationDate(new Date());
                x.setModificationDate(new Date());
                idsIncoming.add(x);
            } else {
                counters.increment("multimatchByRefSeq");

                String info = null;
                for( Transcript tr: transcripts ) {
                    if( info==null ) {
                        info = tr.getAccId()+" (RGD:"+tr.getRgdId()+")";
                    } else {
                        info += "  , "+tr.getAccId()+" (RGD:"+tr.getRgdId()+")";
                    }
                }
                logMultimatch.debug(species+": "+accId+" matches multiple transcripts: "+info);
            }
        }
        in.close();
    }

    void parseRgdFile( List<XdbId> idsIncoming, String taxonId, String species, CounterPool counters ) throws Exception {

        // to avoid inserting duplicate incoming ids
        Set<String> incomingXdbIdsKeys = new HashSet<>();
        for( XdbId xdbId: idsIncoming ) {
            incomingXdbIdsKeys.add(xdbId.getAccId()+"|"+xdbId.getRgdId());
        }

        // file content
        //Tab-separated file with RNAcentral ids, corresponding external ids,
        //NCBI taxon ids, RNA types (according to INSDC classification),
        //        and gene names.
        //URS0000013967	RGD	2325598	10116	pre_miRNA	Mir21
        //URS0000060CC3	RGD	727842	10116	snoRNA	LOC252890

        FileDownloader fd = new FileDownloader();
        fd.setExternalFile(getRgdMappingFile());
        fd.setLocalFile("data/rgd_mapping");
        fd.setUseCompression(true);
        fd.setPrependDateStamp(true);
        String localFile = fd.downloadNew();

        BufferedReader in = Utils.openReader(localFile);
        String line;
        while( (line=in.readLine())!=null ) {
            String[] cols = line.split("[\\t]", -1);
            String rnaCentralId = cols[0];
            String tag = cols[1]; // must be 'RGD'
            String accId = cols[2]; // RGD_ID
            String taxon = cols[3]; // must be '10116'
            String rnaType = cols[4];
            String geneSymbol = cols[5];

            if( !taxonId.equals(taxon) ) {
                continue;
            }
            if( !tag.equals("RGD") ) {
                log.error("*** unexpected db tag: "+tag+";  was expecting: RGD");
                continue;
            }
            counters.increment("linesProcessedForSpecies");

            int rgdId = Integer.parseInt(accId);
            Gene gene = dao.getGeneByRgdId(rgdId);
            if( gene == null ) {

                // todo: consider adding matching by gene symbol
                counters.increment("noMatchByRgdId");
                log.debug("-- no match for " + accId + " gene " + geneSymbol + "  species " + species);

            } else {
                counters.increment("matchByRgdId");

                // do not create duplicate entries
                String key = rnaCentralId+"|"+rgdId;
                if( !incomingXdbIdsKeys.contains(key) ) {
                    XdbId x = new XdbId();
                    x.setAccId(rnaCentralId);
                    x.setSrcPipeline(getPipelineName());
                    x.setRgdId(rgdId);
                    x.setXdbKey(getXdbKeyForRNACentral());
                    x.setCreationDate(new Date());
                    x.setModificationDate(new Date());
                    idsIncoming.add(x);

                    incomingXdbIdsKeys.add(key);
                } else {
                    counters.increment("matchByRgdIdDuplicates");
                }
            }
        }
        in.close();
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setPipelineName(String pipelineName) {
        this.pipelineName = pipelineName;
    }

    public String getPipelineName() {
        return pipelineName;
    }

    public void setXdbKeyForRNACentral(int xdbKeyForRNACentral) {
        this.xdbKeyForRNACentral = xdbKeyForRNACentral;
    }

    public int getXdbKeyForRNACentral() {
        return xdbKeyForRNACentral;
    }

    public void setRefSeqMappingFile(String refSeqMappingFile) {
        this.refSeqMappingFile = refSeqMappingFile;
    }

    public String getRefSeqMappingFile() {
        return refSeqMappingFile;
    }

    public String getRgdMappingFile() {
        return rgdMappingFile;
    }

    public void setRgdMappingFile(String rgdMappingFile) {
        this.rgdMappingFile = rgdMappingFile;
    }
}

