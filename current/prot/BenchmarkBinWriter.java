package prot;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import shared.KillSwitch;
import structures.ByteBuilder;

/**
 * Splits a shred FASTA into per-bin FASTA files using a benchmark manifest (binID, contig, role)
 * emitted by {@link MagQCVectorMaker}'s benchmark mode, so CheckM1/CheckM2 can score the SAME
 * synthetic bins our net scored from vectors. ONE streaming pass over the (large) shred FASTA; a
 * shred appearing in multiple bins is copied to each. Only wanted records accumulate sequence, so
 * the ~99% of shreds no bin uses cost only a name lookup. Crashes loud (KillSwitch) if any manifest
 * contig is absent from the shreds - a stale manifest or wrong FASTA must never silently yield short bins.
 *
 * Usage: java -ea -cp current prot.BenchmarkBinWriter manifest=bench_lean_manifest.tsv \
 *          shreds=foundation_v3/shreds/all_shreds.fa outdir=bins/
 * @author UMP45
 */
public class BenchmarkBinWriter {

	public static void main(String[] args){
		String manifest=null, shreds=null, outdir=null;
		for(final String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("manifest")){manifest=b;}
			else if(a.equals("shreds") || a.equals("in")){shreds=b;}
			else if(a.equals("outdir") || a.equals("out")){outdir=b;}
			else{throw new RuntimeException("Unknown arg: "+arg);}
		}
		if(manifest==null || shreds==null || outdir==null){
			throw new RuntimeException("Usage: manifest=<tsv> shreds=<fasta> outdir=<dir>");
		}
		if(!outdir.endsWith("/")){outdir=outdir+"/";}
		new BenchmarkBinWriter().process(manifest, shreds, outdir);
	}

	private void process(String manifestFile, String shredsFile, String outdir){
		// 1) Load manifest: contig -> bins that want it; open one writer per bin.
		final HashMap<String,ArrayList<String>> want=new HashMap<String,ArrayList<String>>();
		final HashMap<String,ByteStreamWriter> writers=new HashMap<String,ByteStreamWriter>();
		long manifestRows=0;
		{
			final ByteFile bf=ByteFile.makeByteFile(manifestFile, true);
			final LineParser1 lp=new LineParser1((byte)'\t');
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()<2){continue;}
				final String bin=lp.parseString(0);
				final String contig=lp.parseString(1);
				ArrayList<String> l=want.get(contig);
				if(l==null){want.put(contig, l=new ArrayList<String>(1));}
				l.add(bin);
				if(!writers.containsKey(bin)){
					final ByteStreamWriter bsw=new ByteStreamWriter(outdir+bin+".fa", true, false, true);
					bsw.start();
					writers.put(bin, bsw);
				}
				manifestRows++;
			}
			bf.close();
		}
		final int uniqueContigs=want.size();
		assert(manifestRows>0 && uniqueContigs>0) : "empty manifest "+manifestFile;
		System.err.println("manifest: "+manifestRows+" rows, "+uniqueContigs+" unique contigs, "+writers.size()+" bins");

		// 2) Stream shreds; route each WANTED record to its bin(s). Only wanted records buffer sequence.
		final ByteFile bf=ByteFile.makeByteFile(shredsFile, true);
		final ByteBuilder seq=new ByteBuilder(24000);
		final ByteBuilder rec=new ByteBuilder(24200);
		String curName=null; boolean wanted=false;
		long records=0, routed=0; int found=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				if(wanted){routed+=writeRecord(curName, seq, rec, want, writers); found++;}
				curName=contigName(line);
				wanted=want.containsKey(curName);
				if(wanted){seq.clear();}
				records++;
			}else if(wanted){
				seq.append(line);
			}
		}
		if(wanted){routed+=writeRecord(curName, seq, rec, want, writers); found++;}
		bf.close();

		for(final ByteStreamWriter bsw : writers.values()){bsw.poisonAndWait();}
		System.err.println("scanned "+records+" shreds; routed "+routed+" (contig,bin) copies across "+writers.size()+" bins");
		// Crash loud if any manifest contig never appeared in the shreds (wrong FASTA / stale manifest).
		assert(found==uniqueContigs) : KillSwitch.assertDie("BenchmarkBinWriter: "+(uniqueContigs-found)
			+" of "+uniqueContigs+" manifest contigs NOT found in "+shredsFile
			+" - stale manifest or wrong shred FASTA; bins would be short. (found="+found+")");
		assert(routed==manifestRows) : KillSwitch.assertDie("BenchmarkBinWriter: routed "+routed
			+" != manifest rows "+manifestRows+" - a wanted contig was written the wrong number of times");
		System.err.println("done: all "+uniqueContigs+" contigs found; "+routed+" copies == "+manifestRows+" manifest rows.");
	}

	/** Writes ">name\nseq\n" to every bin that wants this contig; returns the number of bins written. */
	private static int writeRecord(String name, ByteBuilder seq, ByteBuilder rec,
			HashMap<String,ArrayList<String>> want, HashMap<String,ByteStreamWriter> writers){
		rec.clear();
		rec.append('>').append(name).nl();
		rec.append(seq.array, 0, seq.length).nl();
		final byte[] recBytes=rec.toBytes();
		final ArrayList<String> bins=want.get(name);
		for(int i=0; i<bins.size(); i++){writers.get(bins.get(i)).print(recBytes);}
		return bins.size();
	}

	/** Contig name = the FULL FASTA header (after '>') with internal spaces canonicalized to '_',
	 *  matching how CacheBuilder names shreds. Most headers are the clean form
	 *  {@code <acc>_tid_<tid>_<coords>} (no spaces -> unchanged); pipe-form headers
	 *  {@code tid|<tid>|<acc> <description>_<coords>} carry a space-delimited description that the
	 *  cache stores with '_' in place of each space. Truncating at the first space (the old bug)
	 *  dropped every descriptive shred, so ~6% of manifest contigs went unmatched. Trailing
	 *  whitespace/CR is stripped so a CRLF FASTA cannot desync the key from the cache. */
	private static String contigName(byte[] header){
		int end=header.length;
		while(end>1 && header[end-1]<=' '){end--;}//strip trailing whitespace/CR
		final byte[] buf=new byte[end-1];
		for(int i=1; i<end; i++){final byte c=header[i]; buf[i-1]=(c==' ' ? (byte)'_' : c);}
		return new String(buf);
	}
}
