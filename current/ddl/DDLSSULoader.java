package ddl;

import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import map.IntObjectMap;
import stream.Read;
import stream.StreamerFactory;

/**
 * Loads SSU (16S/18S) ribosomal sequences from FASTA files and attaches
 * them to DDLRecords by matching TaxIDs.  Headers must contain tid|NNNN
 * or tid_NNNN (parsed by DDLWriter.parseTID).
 *
 * @author Noire, Brian Bushnell
 * @date May 18, 2026
 */
public class DDLSSULoader {

	/** Loads 16S and 18S maps, then attaches to records.  Convenience wrapper. */
	public static void loadAndAttach(ArrayList<DDLRecord> records, String r16sFile, String r18sFile){
		IntObjectMap<byte[]> map16=(r16sFile!=null ? loadSSUMap(r16sFile) : null);
		IntObjectMap<byte[]> map18=(r18sFile!=null ? loadSSUMap(r18sFile) : null);
		attachSSU(records, map16, map18);
	}

	/** Loads SSU from defaults (resources directory). */
	public static void loadAndAttachDefaults(ArrayList<DDLRecord> records){
		IntObjectMap<byte[]>[] maps=loadSSUMapsDefaults();
		attachSSU(records, maps[0], maps[1]);
	}

	/** Loads 16S and 18S maps from default resource paths.
	 * @return {map16S, map18S} — either may be null if resource not found */
	@SuppressWarnings("unchecked")
	public static IntObjectMap<byte[]>[] loadSSUMapsDefaults(){
		String r16s=Data.findPath(DEFAULT_16S);
		String r18s=Data.findPath(DEFAULT_18S);
		if(r16s==null){System.err.println("Warning: 16S file not found: "+DEFAULT_16S);}
		if(r18s==null){System.err.println("Warning: 18S file not found: "+DEFAULT_18S);}
		IntObjectMap<byte[]> map16=(r16s!=null ? loadSSUMap(r16s) : null);
		IntObjectMap<byte[]> map18=(r18s!=null ? loadSSUMap(r18s) : null);
		return new IntObjectMap[]{map16, map18};
	}

	/** Attaches preloaded SSU maps to records by TaxID.
	 * @param records DDL records to attach SSU sequences to
	 * @param map16 TaxID->16S sequence map (null to skip)
	 * @param map18 TaxID->18S sequence map (null to skip) */
	public static void attachSSU(ArrayList<DDLRecord> records,
			IntObjectMap<byte[]> map16, IntObjectMap<byte[]> map18){
		int a16=0, a18=0;
		for(DDLRecord rec : records){
			if(rec.taxID<=0){continue;}
			if(map16!=null && rec.r16S==null){
				byte[] seq=map16.get(rec.taxID);
				if(seq!=null){rec.r16S=seq; a16++;}
			}
			if(map18!=null && rec.r18S==null){
				byte[] seq=map18.get(rec.taxID);
				if(seq!=null){rec.r18S=seq; a18++;}
			}
		}
		System.err.println("Attached "+a16+" 16S and "+a18+" 18S sequences to "+records.size()+" records.");
	}

	private static IntObjectMap<byte[]> loadSSUMap(String path){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		ArrayList<Read> reads=StreamerFactory.getReads(-1, false, ff, null, null, null);
		IntObjectMap<byte[]> map=new IntObjectMap<>(reads.size()*2);
		for(Read r : reads){
			int tid=DDLWriter.parseTID(r.id);
			if(tid>0){map.put(tid, r.bases);}
		}
		return map;
	}

	public static final String DEFAULT_16S="?all_prok_16S_best_taxsorted.fa.gz";
	public static final String DEFAULT_18S="?all_euk_18S_best_taxsorted.fa.gz";
}
