package synth;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import shared.Random;
import shared.Timer;
import stream.Read;

/** Exact truth-CIGAR and disabled-path parity checks for RandomReadsMG.
 * @author Fischl */
public class RandomReadsMGCigarTest {

	public static void main(final String[] args) throws Exception{
		boolean assertions=false;assert(assertions=true);
		if(!assertions){throw new IllegalStateException("Run this regression with -ea.");}
		assert(args.length<=1) : "Optional argument is the BBTools root containing testdata.";
		testLiteralCigars();testMutationParity();testGenericIndelParity();testArguments();
		if(args.length==1){testGeneratedHeaders(args[0]);}
		System.out.println("RANDOMREADSMG_CIGAR_TEST_OK");
	}

	private static void testLiteralCigars(){
		check("ACGT","ACGT",new int[]{0,1,2,3},false,"4=",4);
		check("ACGTAC","ACATAC",new int[]{0,1,-1,3,4,5},false,"2=1I1D3=",6);
		check("ACGTAC","ACATAC",new int[]{0,1,-1,3,4,5},true,"3=1D1I2=",6);
		check("ANCG","ATCG",new int[]{0,1,2,3},false,"1=1M2=",4);
		check("ACG","CG",new int[]{1,2},false,"1D2=",3);
		check("ACGT","AC",new int[]{0,1},false,"2=",2);
		assert(RandomReadsMG.truthStart(100,10,7,false)==100) : "Plus truth starts at the extracted interval start.";
		assert(RandomReadsMG.truthStart(100,10,7,true)==103) : "Minus truth starts at the left edge of its retained reference suffix.";
	}

	private static void check(final String reference,final String observed,final int[] coordinates,
		final boolean reverse,final String expected,final int expectedReference){
		final RandomReadsMG.TruthTrace trace=new RandomReadsMG.TruthTrace(ascii(reference));
		trace.coordinates=coordinates.clone();
		final Read read=new Read(ascii(observed),null,"truth",0);
		final RandomReadsMG.TruthCigar truth=RandomReadsMG.buildTruthCigar(read,trace,reverse);
		assert(Arrays.equals(truth.cigar,ascii(expected))) :
			"Unexpected truth CIGAR: expected="+expected+" actual="+new String(truth.cigar,StandardCharsets.US_ASCII);
		assert(truth.readLength==observed.length() && truth.referenceLength==expectedReference) :
			"Truth consumption mismatch for "+expected;
	}

	private static void testMutationParity(){
		final RandomReadsMG generator=new RandomReadsMG(new String[]{"t=1","tree=f"});
		for(int seed=0;seed<64;seed++){
			final java.util.Random source=new java.util.Random(seed);
			final byte[] input=new byte[257],alphabet=ascii("ACGTN");
			for(int i=0;i<input.length;i++){input[i]=alphabet[source.nextInt(alphabet.length)];}
			final Read legacy=new Read(input.clone(),null,"legacy",0),traced=new Read(input.clone(),null,"traced",0);
			final SeededRandom a=new SeededRandom(seed),b=new SeededRandom(seed);
			final int ca=generator.mutateLongRead(legacy,0.01f,0.02f,0.03f,0.001f,a);
			final RandomReadsMG.TruthTrace trace=new RandomReadsMG.TruthTrace(traced.bases);
			final int cb=generator.mutateLongRead(traced,0.01f,0.02f,0.03f,0.001f,b,trace);
			assert(ca==cb && Arrays.equals(legacy.bases,traced.bases) && a.nextLong()==b.nextLong()) :
				"Truth tracking changed long-read mutation output or RNG state; seed="+seed;
			final RandomReadsMG.TruthCigar truth=RandomReadsMG.buildTruthCigar(traced,trace,false);
			assert(truth.readLength==traced.length()) : "Traced long-read CIGAR lost output bases.";
		}
	}

	private static void testGenericIndelParity(){
		for(int seed=0;seed<64;seed++){
			final byte[] input=new byte[240];Arrays.fill(input,(byte)'A');
			final Read legacy=new Read(input.clone(),null,"legacy",0),traced=new Read(input.clone(),null,"traced",0);
			final SeededRandom a=new SeededRandom(seed),b=new SeededRandom(seed);
			final int ra=RandomReadsMG.addIndels(legacy,0.03f,0.04f,200,30,0,a);
			final RandomReadsMG.TruthTrace trace=new RandomReadsMG.TruthTrace(traced.bases);
			final int rb=RandomReadsMG.addIndels(traced,0.03f,0.04f,200,30,0,b,trace);
			assert(ra==rb && Arrays.equals(legacy.bases,traced.bases) && a.nextLong()==b.nextLong()) :
				"Truth tracking changed generic-indel output or RNG state; seed="+seed;
			final RandomReadsMG.TruthCigar truth=RandomReadsMG.buildTruthCigar(traced,trace,false);
			assert(truth.readLength==200) : "Generic-indel trace did not retain the requested output length.";
		}
	}

	private static void testArguments(){
		new RandomReadsMG(new String[]{"pacbio","addcigar=t","tree=f"});
		expectFailure(new String[]{"addcigar=t","tree=f"});
		expectFailure(new String[]{"pacbio","addcigar=t","circular=t","tree=f"});
	}

	private static void expectFailure(final String[] args){
		boolean failed=false;try{new RandomReadsMG(args);}catch(IllegalArgumentException expected){failed=true;}
		assert(failed) : "Invalid addcigar combination was accepted: "+Arrays.toString(args);
	}

	/** End-to-end three-pass validation against the source FASTA, including minus reads. */
	private static void testGeneratedHeaders(final String root) throws Exception{
		final String ref=root+"/testdata/read_threaded_x_ref.fa";
		final File out=File.createTempFile("randomreadsmg-cigar-", ".fq");
		if(!out.delete()){throw new IllegalStateException("Could not reserve a fresh test output: "+out);}
		try{
			final RandomReadsMG generator=new RandomReadsMG(new String[]{"in="+ref,"out="+out,
				"pacbio","addcigar=t","adderrors=t","srate=0.04","irate=0.03","drate=0.03","hrate=0.02",
				"subrate=0.04","insrate=0.03","delrate=0.04","minlen=40","meanlen=60","maxlen=80",
				"readspercontig=5","seed=23","t=1","tree=f","overwrite=t"});
			generator.process(new Timer());
			final ArrayList<String> references=readFasta(ref);
			int records=0;boolean plus=false,minus=false,sub=false,ins=false,del=false,mixed=false,terminal=false;
			final BufferedReader br=new BufferedReader(new FileReader(out));
			try{
				for(String header;(header=br.readLine())!=null;){
					final String bases=br.readLine();final String plusLine=br.readLine();final String quality=br.readLine();
					assert(bases!=null && "+".equals(plusLine) && quality!=null && bases.length()==quality.length()) :
						"Generated FASTQ must retain four complete lines with matching base/quality lengths.";
					final Matcher hm=HEADER.matcher(header);
					assert(hm.matches()) : "Generated header lacks the exact RandomReadsMG CIGAR schema: "+header;
					final int contig=Integer.parseInt(hm.group(1)),strand=Integer.parseInt(hm.group(2));
					final int start=Integer.parseInt(hm.group(3)),readLength=Integer.parseInt(hm.group(4));
					final int referenceLength=Integer.parseInt(hm.group(5));final String cigar=hm.group(6);
					assert(readLength==bases.length()) : "Header i field must equal FASTQ read length.";
					final String oriented=strand==0 ? bases : reverseComplement(bases);
					final String reference=references.get(contig);
					int readPos=0,refPos=0,lastEnd=0;boolean hasSub=false,hasIndel=false;byte first=0,last=0;
					final Matcher cm=CIGAR.matcher(cigar);
					while(cm.find()){
						assert(cm.start()==lastEnd) : "Unparsed CIGAR text at "+lastEnd+" in "+cigar;
						lastEnd=cm.end();final int length=Integer.parseInt(cm.group(1));final byte op=(byte)cm.group(2).charAt(0);
						if(first==0){first=op;}last=op;
						if(op=='='){
							assert(oriented.regionMatches(readPos,reference,start+refPos,length)) : "= run disagrees with source reference: "+header;
							readPos+=length;refPos+=length;
						}else if(op=='X'){
							for(int i=0;i<length;i++){
								final char a=oriented.charAt(readPos+i),b=reference.charAt(start+refPos+i);
								assert(defined(a) && defined(b) && a!=b) : "X run contains a match or undefined base: "+header;
							}
							readPos+=length;refPos+=length;sub=hasSub=true;
						}else if(op=='M'){
							for(int i=0;i<length;i++){
								assert(!defined(oriented.charAt(readPos+i)) || !defined(reference.charAt(start+refPos+i))) :
									"Ambiguous M must contain an undefined read or reference base: "+header;
							}
							readPos+=length;refPos+=length;
						}else if(op=='I'){readPos+=length;ins=hasIndel=true;}
						else if(op=='D'){refPos+=length;del=hasIndel=true;}
						else{throw new AssertionError("Unsupported truth-CIGAR operation: "+(char)op);}
					}
					assert(lastEnd==cigar.length()) : "Trailing unparsed CIGAR text: "+cigar;
					assert(readPos==readLength && refPos==referenceLength) : "Header/CIGAR consumption mismatch: "+header;
					assert(start>=0 && start+referenceLength<=reference.length()) : "Truth span exceeds source contig: "+header;
					plus|=strand==0;minus|=strand==1;mixed|=hasSub && hasIndel;terminal|=first=='I' || first=='D' || last=='I' || last=='D';records++;
				}
			}finally{br.close();}
			assert(records>0 && plus && minus && sub && ins && del && mixed && terminal) :
				"Integration fixture must cover both strands, S/I/D, mixed passes, and a terminal indel; records="+records;
		}finally{if(out.exists() && !out.delete()){out.deleteOnExit();}}
	}

	private static ArrayList<String> readFasta(final String path) throws Exception{
		final ArrayList<String> list=new ArrayList<String>();StringBuilder current=null;
		final BufferedReader br=new BufferedReader(new FileReader(path));
		try{for(String line;(line=br.readLine())!=null;){
			if(line.startsWith(">")){if(current!=null){list.add(current.toString());}current=new StringBuilder();}
			else if(line.length()>0){assert(current!=null) : "FASTA sequence preceded its header.";current.append(line.toUpperCase());}
		}}finally{br.close();}
		if(current!=null){list.add(current.toString());}assert(!list.isEmpty()) : "Integration reference is empty.";return list;
	}

	private static String reverseComplement(final String s){
		final char[] out=new char[s.length()];
		for(int i=0;i<s.length();i++){
			final char b=s.charAt(s.length()-1-i);out[i]=b=='A'?'T':b=='C'?'G':b=='G'?'C':b=='T'?'A':'N';
		}
		return new String(out);
	}

	private static boolean defined(final char b){return b=='A' || b=='C' || b=='G' || b=='T';}

	private static byte[] ascii(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static final Pattern HEADER=Pattern.compile("^@f_\\d+_c_(\\d+)_s_([01])_p_(\\d+)_i_(\\d+)_r_(\\d+).*_cigar_([^ ]+) 1:$");
	private static final Pattern CIGAR=Pattern.compile("(\\d+)([=XIDM])");

	private static final class SeededRandom implements Random{
		SeededRandom(final long seed){rng=new java.util.Random(seed);}
		@Override public long nextLong(){return rng.nextLong();}
		@Override public long nextLong(final long bound){throw new AssertionError("Unexpected bounded draw.");}
		@Override public double nextGaussian(){return rng.nextGaussian();}
		@Override public void setSeed(final long seed){rng.setSeed(seed);}
		final java.util.Random rng;
	}
}
