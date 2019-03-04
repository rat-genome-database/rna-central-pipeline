package edu.mcw.rgd.rnacentral;

import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.XdbId;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.io.BufferedReader;
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
    private String refSeqMappingFile;

    Logger log = Logger.getLogger("status");

    public static void main(String[] args) throws Exception {

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Main manager = (Main) (bf.getBean("manager"));

        try {
            manager.run();
        }catch (Exception e) {
            manager.log.error(e);
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
        fd.setPrependDateStamp(true);
        String localFile = fd.downloadNew();


        List<Integer> speciesTypeKeys = new ArrayList<>(SpeciesType.getSpeciesTypeKeys());
        Collections.shuffle(speciesTypeKeys);

        for( int speciesTypeKey: speciesTypeKeys ) {
            if( speciesTypeKey!=0 ) {
                int linesProcessed = run(speciesTypeKey, localFile);
                log.info("   "+linesProcessed+" RefSeq lines processed for "+SpeciesType.getCommonName(speciesTypeKey));
            }
        }

        // TODO: download ensembl file and process it as well

        msg = "===    time elapsed: "+ Utils.formatElapsedTime(startTime, System.currentTimeMillis());
        log.info(msg);

        log.info("");
    }

    public int run(int speciesTypeKey, String localFile) throws Exception {

        // get taxon id for given species
        final String taxonId = Integer.toString(SpeciesType.getTaxonomicId(speciesTypeKey));
        int linesProcessedForSpecies = 0;

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
            linesProcessedForSpecies++;

            System.out.println(line);
        }






        String species = SpeciesType.getCommonName(speciesTypeKey);
        String msg = "START: " + getPipelineName() + " ID generation starting for " + species;
        log.info(msg);

        // QC
        log.debug("QC: get "+getPipelineName()+" Ids in RGD for "+species);
        List<XdbId> idsInRgd = dao.getHumanProteomeMapIds(speciesTypeKey, getPipelineName());
        log.debug("QC: get incoming "+getPipelineName()+" Ids for "+species);
        List<XdbId> idsIncoming = getIncomingIds(speciesTypeKey);

        // determine to-be-inserted Human Proteome Map ids
        log.debug("QC: determine to-be-inserted "+getPipelineName()+" Ids");
        List<XdbId> idsToBeInserted = new ArrayList<XdbId>(idsIncoming);
        idsToBeInserted.removeAll(idsInRgd);

        // determine matching Human Proteome Map ids
        log.debug("QC: determine matching "+getPipelineName()+" Ids");
        List<XdbId> idsMatching = new ArrayList<XdbId>(idsIncoming);
        idsMatching.retainAll(idsInRgd);

        // determine to-be-deleted Human Proteome Map ids
        log.debug("QC: determine to-be-deleted "+getPipelineName()+" Ids");
        idsInRgd.removeAll(idsIncoming);
        List<XdbId> idsToBeDeleted = idsInRgd;


        // loading
        if( !idsToBeInserted.isEmpty() ) {
            msg = "  inserting "+getPipelineName()+" ids for "+species+": "+idsToBeInserted.size();
            log.info(msg);
            dao.insertXdbs(idsToBeInserted);
        }

        if( !idsToBeDeleted.isEmpty() ) {
            msg = "  deleting "+getPipelineName()+" ids for "+species+": "+idsToBeDeleted.size();
            log.info(msg);
            dao.deleteXdbIds(idsToBeDeleted);
        }

        if( !idsMatching.isEmpty() ) {
            msg = "  matching "+getPipelineName()+" ids for "+species+": "+idsMatching.size();
            log.info(msg);
            dao.updateModificationDate(idsMatching);
        }

        msg = "END: "+getPipelineName() + " ID generation complete for " + species;
        log.info(msg);

        return linesProcessedForSpecies;
    }

    List<XdbId> getIncomingIds(int speciesTypeKey) throws Exception {

        List<Gene> genes = dao.getActiveGenes(speciesTypeKey);
        List<XdbId> incomingIds = new ArrayList<XdbId>(genes.size());
        for (Gene g: genes) {
            XdbId x = new XdbId();
            x.setAccId(g.getSymbol());
            x.setSrcPipeline(getPipelineName());
            x.setRgdId(g.getRgdId());
            x.setXdbKey(56);
            x.setCreationDate(new Date());
            x.setModificationDate(new Date());
            incomingIds.add(x);
        }
        return incomingIds;
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

    public void setRefSeqMappingFile(String refSeqMappingFile) {
        this.refSeqMappingFile = refSeqMappingFile;
    }

    public String getRefSeqMappingFile() {
        return refSeqMappingFile;
    }
}

