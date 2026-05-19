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

	/** Loads 16S and/or 18S sequences and attaches them to records by TaxID.
	 * @param records DDL records to attach SSU sequences to
	 * @param r16sFile Path to 16S FASTA (null to skip)
	 * @param r18sFile Path to 18S FASTA (null to skip) */
	public static void loadAndAttach(ArrayList<DDLRecord> records, String r16sFile, String r18sFile){
		int a16=0, a18=0;
		if(r16sFile!=null){
			IntObjectMap<byte[]> map16=loadSSUMap(r16sFile);
			for(DDLRecord rec : records){
				if(rec.taxID>0 && rec.r16S==null){
					byte[] seq=map16.get(rec.taxID);
					if(seq!=null){rec.r16S=seq; a16++;}
				}
			}
		}
		if(r18sFile!=null){
			IntObjectMap<byte[]> map18=loadSSUMap(r18sFile);
			for(DDLRecord rec : records){
				if(rec.taxID>0 && rec.r18S==null){
					byte[] seq=map18.get(rec.taxID);
					if(seq!=null){rec.r18S=seq; a18++;}
				}
			}
		}
		System.err.println("Attached "+a16+" 16S and "+a18+" 18S sequences to "+records.size()+" records.");
	}

	/** Loads SSU from defaults (resources directory). */
	public static void loadAndAttachDefaults(ArrayList<DDLRecord> records){
		String r16s=Data.findPath(DEFAULT_16S);
		String r18s=Data.findPath(DEFAULT_18S);
		if(r16s==null){System.err.println("Warning: 16S file not found: "+DEFAULT_16S);}
		if(r18s==null){System.err.println("Warning: 18S file not found: "+DEFAULT_18S);}
		loadAndAttach(records, r16s, r18s);
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
