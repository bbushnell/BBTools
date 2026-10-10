package prot;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import parse.Parser;
import shared.Shared;
import fileIO.ReadWrite;

/** Isolated loader timing and exact native-versus-text model comparison. @author Collei */
public final class HbmDenseTextProbe {

	public static void main(final String[] args) throws Exception{
		String text=null, binary=null, reference=null, provenance=null, pin=null, mode="text", output=null;
		int threads=1, repeats=3, minCount=1;
		boolean phases=false;
		for(final String arg : Parser.parseConfig(args)){
			if(arg.equalsIgnoreCase("selftest=t")){HbmDenseTextLoaderTest.main(new String[0]); return;}
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in")){text=value;}
			else if(key.equals("binary")){binary=value;}
			else if(key.equals("ref")){reference=value;}
			else if(key.equals("provenance")){provenance=value;}
			else if(key.equals("sourcepin")){pin=value;}
			else if(key.equals("mode")){mode=value;}
			else if(key.equals("out")){output=value;}
			else if(key.equals("t")){threads=Integer.parseInt(value);}
			else if(key.equals("repeats")){repeats=Integer.parseInt(value);}
			else if(key.equals("mincount")){minCount=Integer.parseInt(value);}
			else if(key.equals("phases")){phases=parse.Parse.parseBoolean(value);}
			else if(key.equals("blockinput")){HbmCompactTextReader.BLOCK_INPUT=parse.Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		if(mode.equals("decompress")){
			if(text==null || threads<1 || repeats<1){throw new IllegalArgumentException("decompress requires in= t= repeats=");}
			Shared.setThreads(threads); drain(text, threads, repeats); return;
		}
		if(mode.equals("pack")){
			if(minCount!=1){throw new IllegalArgumentException("pack must preserve all stored counts");}
			if(binary==null || reference==null || provenance==null || output==null){throw new IllegalArgumentException("pack requires binary= ref= provenance= out=fresh.hbmt");}
			HbmDenseTextPacker.pack(Paths.get(binary), reference, provenance, Paths.get(output));
			System.out.println("VERIFIED_DENSE_TEXT_PACK_PASS"); return;
		}
		if(text==null || binary==null || reference==null || provenance==null || pin==null ||
				threads<1 || repeats<1 || minCount<1 || !(mode.equals("text") || mode.equals("verified") || mode.equals("binary") || mode.equals("compare"))){
			throw new IllegalArgumentException("in=text binary=native-or-A48 ref=fasta provenance=tsv sourcepin=sha80 mode=text|binary|compare t=1 repeats=3 required");
		}
		Shared.setThreads(threads);
		if(phases && (!mode.equals("binary") || !HbmCompactTextReader.matches(Paths.get(binary)) || minCount!=1)){
			throw new IllegalArgumentException("phases=t requires compact binary= with mincount=1");
		}
		final List<ProteinSequence> sequences=ProteinSearch.readFasta(reference);
		final ArrayList<String> roster=new ArrayList<String>(sequences.size());
		final HashMap<String, byte[]> consensus=new HashMap<String, byte[]>();
		for(final ProteinSequence seq : sequences){
			roster.add(seq.id);
			if(consensus.put(seq.id, seq.enc)!=null){throw new IllegalArgumentException("Duplicate consensus "+seq.id);}
		}
		final byte[][] trusted=HbmBundleLoader.loadSemanticProvenance(provenance);
		if(minCount!=1 && !HbmDenseTextBundle.matches(Paths.get(text))){throw new IllegalArgumentException("Filtered text comparisons require verified v2");}
		if(mode.equals("compare")){
			final HbmBundleLoader.Loaded nativeModels=HbmBundleLoader.load(Paths.get(binary), roster, consensus::get, trusted, minCount);
			final HbmBundleLoader.Loaded textModels=HbmDenseTextBundle.matches(Paths.get(text)) ?
				HbmBundleLoader.load(Paths.get(text), roster, consensus::get, trusted, minCount) :
				HbmDenseTextLoader.load(text, roster, consensus::get, pin, threads);
			textModels.assertStructuralEquivalent(nativeModels);
			System.out.println("EXACT_NATIVE_TEXT_MODELS_PASS\tfamilies="+nativeModels.familyCount()+"\tthreads="+threads+"\tminCount="+minCount);
			return;
		}
		System.out.println("mode\tthreads\trepeat\tfamilies\tload_seconds"+(phases ?
			"\tread_headers_checksum_s\tparse_worker_sum_s\treconstruct_worker_sum_s\tinput_lines_or_framing_s\theader_parse_s\tchecksum_s\tworker_split_sum_s\tblock_stream_read_s" : ""));
		for(int repeat=1; repeat<=repeats; repeat++){
			// GC is outside the loader clock and makes each load reclaimable. This is
			// an isolated loader measurement, not whole-ProkCC startup performance.
			System.gc();
			final HbmCompactTextReader.Timings clocks=phases ? new HbmCompactTextReader.Timings() : null;
			final long start=System.nanoTime();
			final HbmBundleLoader.Loaded models=phases ? HbmCompactTextReader.read(binary, threads, clocks).models() : mode.equals("text") ?
				HbmDenseTextLoader.load(text, roster, consensus::get, pin, threads) :
				HbmBundleLoader.load(Paths.get(mode.equals("verified") ? text : binary), roster, consensus::get, trusted, minCount);
			final double seconds=(System.nanoTime()-start)*1e-9;
			System.out.printf(Locale.ROOT, "%s\t%d\t%d\t%d\t%.6f", mode, threads, repeat, models.familyCount(), seconds);
			if(phases){System.out.printf(Locale.ROOT, "\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f",
				clocks.readNanos*1e-9, clocks.parseNanos*1e-9, clocks.reconstructNanos*1e-9,
				clocks.lineNanos*1e-9, clocks.headerNanos*1e-9, clocks.hashNanos*1e-9, clocks.splitNanos*1e-9, clocks.rawIoNanos*1e-9);}
			System.out.println();
		}
	}

	/** Isolated native decompression plus file I/O; excludes lines, checksums, parsing and graphs. */
	private static void drain(String path, int threads, int repeats) throws Exception{
		final byte[] buffer=new byte[65536]; long expected=-1;
		System.out.println("mode\tthreads\trepeat\tdecompressed_bytes\tload_seconds");
		for(int repeat=1; repeat<=repeats; repeat++){
			System.gc(); long bytes=0; final long start=System.nanoTime();
			final java.io.InputStream input=ReadWrite.getInputStream(path, true, false);
			try{
				for(int n=input.read(buffer); n>=0; n=input.read(buffer)){bytes+=n;}
			}finally{if(ReadWrite.finishReading(input, path, false)){throw new java.io.IOException("Decompression input failed: "+path);}}
			final double seconds=(System.nanoTime()-start)*1e-9;
			if(bytes<1 || (expected>=0 && bytes!=expected)){throw new AssertionError("Repeated decompression changed byte count: "+bytes+" vs "+expected);}
			expected=bytes;
			System.out.printf(Locale.ROOT, "decompress\t%d\t%d\t%d\t%.6f%n", threads, repeat, bytes, seconds);
		}
	}
	private HbmDenseTextProbe(){}
}
