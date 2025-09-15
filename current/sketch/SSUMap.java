package sketch;

import java.io.PrintStream;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;
import tax.GiToTaxid;
import tax.TaxTree;

public class SSUMap {
	
	public static synchronized void load(PrintStream outstream){
		if(r16SFile!=null && r16SMap==null){
			r16SMap=load(r16SFile, TaxTree.default16SFile(), outstream);
		}
		if(r18SFile!=null && r18SMap==null){
			r18SMap=load(r18SFile, TaxTree.default18SFile(), outstream);
		}
	}
	
	private static synchronized HashMap<Integer, byte[]> load(String ssuFile, String defaultFile, PrintStream outstream){
		HashMap<Integer, byte[]> map=null;
		if(ssuFile!=null){
			final boolean oldAminoIn=Shared.AMINO_IN;
			final boolean oldInterleaved=FASTQ.FORCE_INTERLEAVED;
			Shared.AMINO_IN=false;//SSUs are nucleotide, which can cause a crash, esp. with IUPAC symbols
			FASTQ.FORCE_INTERLEAVED=false;

			String fname=ssuFile;
			if("auto".equalsIgnoreCase(fname)){fname=defaultFile;}
			final FileFormat ffssu=FileFormat.testInput(fname, FileFormat.FA, null, true, false);
			map=loadSSU(ffssu, outstream);

			Shared.AMINO_IN=oldAminoIn;
			FASTQ.FORCE_INTERLEAVED=oldInterleaved;
		}
		return map;
	}
	
	private static HashMap<Integer, byte[]> loadSSU(FileFormat ff, PrintStream outstream){
		ConcurrentReadInputStream cris=makeCris(ff, outstream);
		HashMap<Integer, byte[]> map=new HashMap<Integer, byte[]>(1000000);
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();

		//Check to ensure pairing is as expected
		if(ln!=null && !ln.isEmpty()){
//			if(verbose){outstream.println("Fetched "+ln.size()+" reads.");}
			Read r=ln.get(0);
			assert(ff.samOrBam() || (r.mate!=null)==cris.paired());
		}

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
			if(verbose){outstream.println("Fetched "+ln.size()+" reads.");}
			for(Read r : ln){
				final int tid=GiToTaxid.getID(r.id);
				if(tid>=0 && r.length()>1000){
					byte[] old=map.get(tid);
					if(old==null || old.length<r.length()){map.put(tid, r.bases);}
				}
			}
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());

			//Fetch a new list
			ln=cris.nextList();
		}

		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		errorState|=ReadWrite.closeStream(cris);
		
		return map;
	}
	
	private static ConcurrentReadInputStream makeCris(FileFormat ff, PrintStream outstream){
		if(verbose){outstream.println("makeCris");}
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Loading "+ff.name());}
		boolean paired=cris.paired();
		assert(!paired);
//		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		return cris;
	}

	public static boolean hasMap(){return r16SMap!=null || r18SMap!=null;}
	public static int r16SCount(){return r16SMap==null ? 0 : r16SMap.size();}
	public static int r18SCount(){return r18SMap==null ? 0 : r18SMap.size();}
	
	public static String r16SFile=null;
	public static String r18SFile=null;
	public static HashMap<Integer, byte[]> r16SMap=null;
	public static HashMap<Integer, byte[]> r18SMap=null;
	static boolean verbose=false;
	static boolean errorState=false;
	
}
