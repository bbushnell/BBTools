package test.cardinality;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongList;

/** Compares native/shared-iterator emissions with an independent reset-counter
 * oracle and the preserved pre-fix loop. The legacy loop is test evidence only.
 * Optional real FASTA input is read by the native Streamer without transformation.
 * @author Ganyu
 */
public final class HashSmallAmbiguityTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HashSmallAmbiguityTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			boolean fixed=true;
			boolean requireDefined=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("expect")){
					require("fixed".equals(b) || "legacy".equals(b), "expect=fixed|legacy identifies the actual build under test");
					fixed=b.equals("fixed");
				}else if(a.equals("requiredefined")){requireDefined=Parse.parseBoolean(b);}
				else{require(parser.parse(arg, a, b), "Unknown ambiguity test option: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && parser.out1!=null && parser.out2==null && parser.in2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh.tsv, optional in=real.fa, -ea, no sampling, paired input or overwrite");
			RrnaResourceIO.fresh(parser.out1, "");
			final int exponent=DynamicDemiLog.exponentBits();
			DynamicDemiLog.setExponent(4);
			try(Output out=new Output(parser.out1)){
				out.writer.println("source\tk\tbases\tambiguous_bases\tlegacy_emissions\tvalid_emissions\tremoved_spanning_emissions\tactual_emissions\texpectation");
				final HashSmallAmbiguityTest test=new HashSmallAmbiguityTest(out, fixed);
				test.synthetic();
				if(parser.in1!=null){test.real(parser.in1, requireDefined);}
				pp.outstream.println("HASH_SMALL_AMBIGUITY_TEST_PASS\tcases="+test.cases+"\texpect="+(fixed ? "fixed" : "legacy"));
			}finally{DynamicDemiLog.setExponent(exponent);}
		}catch(Throwable failure){
			// Streamer/writer workers must never keep a failed main-thread test alive.
			failure.printStackTrace(pp.outstream);
			KillSwitch.kill("HashSmall ambiguity test failed");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private HashSmallAmbiguityTest(final Output out_, final boolean fixed_){out=out_; fixed=fixed_;}

	/** Every supported short k, each single-N position, consecutive Ns, lowercase,
	 * all four homopolymers, empty/short prefixes and deterministic mixed sequences. */
	private void synthetic(){
		final Random random=new Random(20261113);
		for(int k=1; k<=31; k++){
			verify("empty", new byte[0], k);
			for(final char base:"ACTGactg".toCharArray()){
				final byte[] bases=new byte[2*k+3];
				Arrays.fill(bases, (byte)base);
				verify("homopolymer_"+base, bases, k);
				verify("short_"+base, Arrays.copyOf(bases, k-1), k);
				for(int pos=0; pos<bases.length; pos++){
					bases[pos]='N';
					verify("single_N_"+base+"_"+pos, bases, k);
					bases[pos]=(byte)base;
				}
			}
			final byte[] boundary=new byte[2*k-1];
			Arrays.fill(boundary, (byte)'A');
			boundary[k-1]='N';
			verify("literal_A_N_A", boundary, k);
			final byte[] allN=new byte[4*k+5];
			Arrays.fill(allN, (byte)'N');
			verify("all_N", allN, k);
			for(int repeat=0; repeat<8; repeat++){
				final byte[] mixed=new byte[512];
				for(int i=0; i<mixed.length; i++){mixed[i]=(byte)"ACGTNNacgtn".charAt(random.nextInt(11));}
				verify("mixed_"+repeat, mixed, k);
			}
		}
	}

	/** Actual sequence data, no masking or concatenation. Require-defined mode
	 * rejects a supposedly N-free genome if the actual input contains an N. */
	private void real(final String path, final boolean requireDefined){
		assert(path!=null) : "Real-file parity requires a named immutable input";
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		final Streamer input=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		input.start();
		long reads=0;
		for(ListNum<Read> list=input.nextList(); list!=null; list=input.nextList()){
			if(list.list==null || list.list.isEmpty()){input.returnList(list); break;}
			for(final Read read:list.list){
				require(read!=null && read.bases!=null && read.mate==null, "Real parity consumes unpaired FASTA records");
				if(requireDefined){require(ambiguous(read.bases)==0, "requiredefined=t input contains ambiguous bases: "+read.id);}
				for(final int k:new int[]{4,10,31}){verify("real_"+(reads+1), read.bases, k);}
				reads++;
			}
			input.returnList(list);
		}
		require(!ReadWrite.closeStreams(input), "Real genome stream did not close cleanly: "+path);
		require(reads>0, "Real genome parity cannot pass on empty input: "+path);
	}

	/** Ordered word equality is stronger than equality of compressed sketch keys.
	 * Filtering old EMISSIONS preserves a valid occurrence of a word even if an
	 * invalid N-spanning occurrence happens to encode the same word elsewhere. */
	private void verify(final String source, final byte[] bases, final int k){
		assert(k>=1 && k<=31) : "The short-word dispatch domain is k1..31";
		final int ns=ambiguous(bases);
		final LongList oracle=oracle(bases, k);
		final Legacy old=legacy(bases, k);
		equal(old.valid, oracle, source, k, "Removing N-spanning old emissions must preserve every defined word");
		final Recorder nativeWords=new Recorder(k, false);
		nativeWords.hashSmall(bases, null);
		equal(nativeWords.words, fixed ? oracle : old.all, source, k, "Actual shared loop must match the declared build expectation");
		if(ns==0){equal(nativeWords.words, old.all, source, k, "N-free words must remain byte-identical to the old loop");}
		final Recorder quality=new Recorder(k, true);
		final byte[] qualities=new byte[bases.length];
		Arrays.fill(qualities, (byte)40);
		quality.hashSmall(bases, qualities);
		equal(quality.words, oracle, source, k, "Quality-aware branch must match the independent counter-reset oracle");
		final Recorder shared=new Recorder(k, false);
		CardinalityTracker.forEachKmer(bases, k, k, true, false, shared);
		equal(shared.words, nativeWords.words, source, k, "The shared iterator must preserve the native loop's exact ordered emissions");
		row.clear().append(source).tab().append(k).tab().append(bases.length).tab().append(ns).tab()
			.append(old.all.size).tab().append(oracle.size).tab().append(old.all.size-old.valid.size).tab()
			.append(nativeWords.words.size).tab().append(fixed ? "fixed" : "legacy").nl();
		out.writer.print(row);
		cases++;
	}

	/** Simple independent reference: reset on N, count defined bases explicitly. */
	private static LongList oracle(final byte[] bases, final int k){
		final LongList words=new LongList();
		final long mask=(1L<<(2*k))-1;
		long forward=0;
		long reverse=0;
		int valid=0;
		for(final byte base:bases){
			if(base=='N' || base=='n'){
				valid=0;
				forward=0;
				reverse=0;
				continue;
			}
			final int code=(base>>1)&3;
			forward=((forward<<2)|code)&mask;
			reverse=(reverse>>>2)|((long)(code^2)<<(2*k-2));
			valid++;
			if(valid>=k){words.add(Math.max(forward, reverse));}
		}
		return words;
	}

	/** Literal original arithmetic from CardinalityTracker at 3472c0d8:601,657–675.
	 * The extra reset counter classifies emitted occurrences; it does not alter them. */
	private static Legacy legacy(final byte[] bases, final int k){
		final Legacy result=new Legacy();
		final int lenmask=(-1)>>>k;
		final int shift2=2*k-2;
		final long mask=(1L<<(2*k))-1;
		long forward=0;
		long reverse=0;
		int len=-1;
		int valid=0;
		for(final byte b:bases){
			final long x=((b>>1)&3)|((b&8)<<28);
			len|=x;
			len>>>=1;
			forward=((forward<<2)|x)&mask;
			reverse=((reverse>>>2)|((x^2)<<shift2))&mask;
			valid=(b=='N' || b=='n') ? 0 : valid+1;
			if(len<=lenmask){
				final long word=Math.max(forward, reverse);
				result.all.add(word);
				if(valid>=k){result.valid.add(word);}
			}
		}
		return result;
	}

	private static int ambiguous(final byte[] bases){
		int n=0;
		for(final byte b:bases){
			if(b=='N' || b=='n'){n++;}
			else if(b!='A' && b!='C' && b!='G' && b!='T' && b!='a' && b!='c' && b!='g' && b!='t'){
				throw new IllegalArgumentException("This N-boundary fixture requires ACGTN; unsupported byte="+b);
			}
		}
		return n;
	}

	private static void equal(final LongList a, final LongList b, final String source, final int k, final String why){
		require(a.size==b.size, why+"; source="+source+", k="+k+", counts="+a.size+"/"+b.size);
		for(int i=0; i<a.size; i++){
			if(a.array[i]!=b.array[i]){
				throw new IllegalArgumentException(why+"; source="+source+", k="+k+", emission="+i);
			}
		}
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Legacy{
		final LongList all=new LongList();
		final LongList valid=new LongList();
	}
	private static final class Recorder extends CardinalityTracker implements CardinalityTracker.KmerConsumer{
		Recorder(final int k, final boolean quality){super(32, k, 12345, quality ? 0.1f : 0f);}
		@Override
		public void accept(final long word, final int length){
			assert(length==k) : "The shared ambiguity fixture must receive only its declared complete word length";
			add(word);
		}
		@Override public void hashAndStore(final long word){words.add(word);}
		@Override public CardinalityTracker copy(){throw new UnsupportedOperationException("Recorder does not estimate");}
		@Override public long cardinality(){throw new UnsupportedOperationException("Recorder does not estimate");}
		@Override public void add(final CardinalityTracker other){throw new UnsupportedOperationException("Recorder does not merge");}
		@Override public float[] compensationFactorLogBucketsArray(){throw new UnsupportedOperationException("Recorder does not calibrate");}
		final LongList words=new LongList();
	}
	private final Output out;
	private final boolean fixed;
	private final ByteBuilder row=new ByteBuilder();
	private long cases;
}
