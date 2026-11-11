package test.assemble;

import java.io.File;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Generates an exhaustive small DNA-word input and independently counts compressed output.
 * Each generated record contains exactly one kmer. The checker expects every raw word,
 * or every reverse-complement equivalence class, exactly once across output records.
 * It does not invoke the compressor or use its table/ownership implementation.
 * @author Brian Bushnell, Shinobu
 */
public final class KmerCompressorCensusTest{

	/** Selects mode=generate out=NEW_FASTA or mode=check in=FASTA, with k=1..12.
	 * rcomp controls the checker's equivalence relation; generation always emits all words. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, KmerCompressorCensusTest.class, false);
		try{
			String mode=null, in=null, out=null;
			int k=8;
			boolean rcomp=true;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected name=value: "+arg);}
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair[1];
				if(key.equals("mode")){mode=value;}
				else if(key.equals("in")){in=value;}
				else if(key.equals("out")){out=value;}
				else if(key.equals("k")){k=Integer.parseInt(value);}
				else if(key.equals("rcomp")){rcomp=Parse.parseBoolean(value);}
				else{throw new IllegalArgumentException("Unknown census option: "+arg);}
			}
			if(k<1 || k>12){throw new IllegalArgumentException("Census k must be in [1,12]: "+k);}
			if("generate".equals(mode) && in==null && out!=null){generate(out, k);}
			else if("check".equals(mode) && in!=null && out==null){check(in, k, rcomp);}
			else{throw new IllegalArgumentException("Use mode=generate out=NEW_FASTA or mode=check in=FASTA.");}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private KmerCompressorCensusTest(){}

	/** Emits all 4^k words in lexical numeric order, one per FASTA record.
	 * Refuses an existing path; a failed write may leave a partial diagnostic fixture. */
	private static void generate(final String path, final int k){
		assert(k>=1 && k<=12) : "The exhaustive generator must fit its bounded integer word domain.";
		if(new File(path).exists()){throw new IllegalArgumentException("Output already exists: "+path);}
		final int words=1<<(2*k);
		final ByteBuilder buffer=new ByteBuilder(65536);
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		final boolean error;
		try{
			for(int word=0; word<words; word++){
				buffer.append('>').append(word).append('\n');
				for(int shift=2*(k-1); shift>=0; shift-=2){buffer.append(BASES[(word>>>shift)&3]);}
				buffer.append('\n');
				if(buffer.length()>=65536){writer.print(buffer); buffer.clear();}
			}
			if(buffer.length()>0){writer.print(buffer);}
		}finally{error=writer.poisonAndWait();}
		if(error){throw new IllegalStateException("Failed to write exhaustive census input: "+path);}
		System.out.println("KCOMPRESS_CENSUS_GENERATED k="+k+" records="+words);
	}

	/** Counts output windows with an independent rolling two-bit encoding.
	 * Headers and N separators break windows; wrapped sequence lines do not.
	 * Saves the complete count summary before rejecting missing or repeated words. */
	private static void check(final String path, final int k, final boolean rcomp){
		assert(k>=1 && k<=12) : "The census bitmap must fit the bounded exhaustive domain.";
		final int words=1<<(2*k), mask=words-1, shift=2*(k-1);
		final byte[] seen=new byte[words];
		int forward=0, reverse=0, run=0, distinct=0;
		long total=0, duplicates=0;
		boolean header=false;
		final ByteFile reader=ByteFile.makeByteFile(path, false);
		final boolean error;
		try{
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){header=true; forward=reverse=run=0; continue;}
				if(!header){throw new IllegalArgumentException("Sequence before FASTA header: "+path);}
				for(byte base : line){
					final int lower=base|32;
					if(lower=='n'){forward=reverse=run=0; continue;}
					final int code=lower=='a' ? 0 : lower=='c' ? 1 : lower=='g' ? 2 : lower=='t' ? 3 : -1;
					if(code<0){throw new IllegalArgumentException("Unexpected FASTA sequence byte: "+base);}
					forward=((forward<<2)|code)&mask;
					reverse=(reverse>>>2)|((3-code)<<shift);
					if(run<k){run++;}
					if(run==k){
						final int word=rcomp ? Math.max(forward, reverse) : forward;
						total++;
						if(seen[word]==0){seen[word]=1; distinct++;}
						else{duplicates++;}
					}
				}
			}
		}finally{error=reader.close();}
		if(error){throw new IllegalStateException("Failed to read census output: "+path);}
		// Odd words cannot equal their DNA reverse complement; even words have 2^k palindromes.
		final int palindromes=(k&1)==0 ? 1<<k : 0;
		final int expected=rcomp ? (words+palindromes)/2 : words;
		final int missing=expected-distinct;
		assert(missing>=0 && total==distinct+duplicates) :
				"Canonical encoding must remain within the exhaustive domain and count every window once.";
		System.out.println("k\trcomp\texpected\ttotal\tdistinct\tmissing\tduplicates");
		System.out.println(k+"\t"+rcomp+"\t"+expected+"\t"+total+"\t"+distinct+"\t"+missing+"\t"+duplicates);
		if(missing!=0 || duplicates!=0){throw new IllegalStateException("Kmer census found missing or repeated words.");}
		System.out.println("KCOMPRESS_CENSUS_PASS");
	}

	/** Independent DNA encoding used only for exhaustive fixture generation. */
	private static final byte[] BASES={'A', 'C', 'G', 'T'};
}
