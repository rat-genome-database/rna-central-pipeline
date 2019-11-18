package edu.mcw.rgd.rnacentral;

import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.datamodel.XdbId;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * @author mtutaj
 * @since 9/3/13
 * <p>
 * wrapper to handle all DAO code
 */
public class DAO {

    XdbIdDAO xdao = new XdbIdDAO();
    TranscriptDAO tdao = new TranscriptDAO();

    Logger logInserted = Logger.getLogger("insertedIds");
    Logger logDeleted = Logger.getLogger("deletedIds");

    public String getConnectionInfo() {
        return xdao.getConnectionInfo();
    }

    public List<Transcript> getTranscriptsByAccId(String accId) throws Exception {
        return tdao.getTranscriptsByAccId(accId);
    }

    public List<XdbId> getRNACentralIds(int speciesTypeKey, String srcPipeline, int xdbKey) throws Exception {

        XdbId filter = new XdbId();
        filter.setXdbKey(xdbKey);
        filter.setSrcPipeline(srcPipeline);
        return xdao.getXdbIds(filter, speciesTypeKey);
    }

    /**
     * insert a bunch of XdbIds; duplicate entries are not inserted (with same RGD_ID,XDB_KEY,ACC_ID,SRC_PIPELINE)
     * @param xdbs list of XdbIds objects to be inserted
     * @return number of actually inserted rows
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int insertXdbs(List<XdbId> xdbs) throws Exception {

        for( XdbId xdbId: xdbs ) {
            logInserted.debug(xdbId.dump("|"));
        }

        return xdao.insertXdbs(xdbs);
    }

    /**
     * delete a list external ids (RGD_ACC_XDB rows);
     * if ACC_XDB_KEY is provided, it is used to delete the row;
     * else ACC_ID, RGD_ID, XDB_KEY and SRC_PIPELINE are used to locate and delete every row
     *
     * @param xdbIds list of external ids to be deleted
     * @return nr of rows deleted
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int deleteXdbIds( List<XdbId> xdbIds ) throws Exception {

        for( XdbId xdbId: xdbIds ) {
            logDeleted.debug(xdbId.dump("|"));
        }

        return xdao.deleteXdbIds(xdbIds);
    }

    public int updateModificationDate(List<XdbId> xdbIds) throws Exception {

        List<Integer> xdbKeys = new ArrayList<Integer>(xdbIds.size());
        for( XdbId xdbId: xdbIds ) {
            xdbKeys.add(xdbId.getKey());
        }
        return xdao.updateModificationDate(xdbKeys);
    }
}
