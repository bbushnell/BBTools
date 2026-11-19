package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.Shared;
import structures.ByteBuilder;
import structures.LongList;

/** Records words emitted by the actual shared hashSmall loop at an N boundary.
 * The recording subclass overrides storage only, never word enumeration.
 * No genome reads, shared-class changes or correction of legacy semantics.
 * @author Ganyu
 */
public final class NcrnaNativeAmbiguityTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaNativeAmbiguityTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				require(parser.parse(arg, a, b), "Unknown native ambiguity fixture option: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && parser.out1!=null && parser.out2==null && parser.in1==null && parser.in2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_evidence.tsv with assertions, no input, sampling or overwrite");
			RrnaResourceIO.fresh(parser.out1, "");
			final int previous=DynamicDemiLog.exponentBits();
			DynamicDemiLog.setExponent(4);
			try(Output out=new Output(parser.out1)){
				out.writer.println("case\troute\tk\tsequence\texpected_valid_words\temitted_words\tordinal\tencoded_word\tdecoded_canonical_word");
				final Recorder shortControl=measure(out, "short_control", "AAA", 4, false, 0);
				final Recorder completeControl=measure(out, "complete_control", "AAAA", 4, false, 1);
				require(shortControl.words.size==0 && completeControl.words.size==1,
					"Actual shared enumeration must distinguish fewer than k defined bases from one complete word");
				int emitted=0;
				for(final int k:new int[]{4, 10}){
					final char[] bases=new char[2*k-1];
					Arrays.fill(bases, 'A');
					bases[k-1]='N';
					final String sequence=new String(bases);
					final Recorder nativeWords=measure(out, "boundary", sequence, k, false, 0);
					final Recorder qualityWords=measure(out, "boundary_quality_control", sequence, k, true, 0);
					require(qualityWords.words.size==0, "The independent quality-aware native branch must reject every N-spanning word");
					final Recorder sharedWords=new Recorder(k, false);
					CardinalityTracker.forEachKmer(sequence.getBytes(StandardCharsets.US_ASCII), k, k, true, false, sharedWords);
					write(out, "boundary_shared_iterator", "shared_iterator_canonical", sequence, k, 0, sharedWords);
					require(Arrays.equals(nativeWords.words.toArray(), sharedWords.words.toArray()),
						"The shared canonical iterator must retain the actual native loop's exact emission order, including legacy ambiguity behavior");
					emitted+=nativeWords.words.size;
				}
				pp.outstream.println(emitted>0 ? "NATIVE_N_BOUNDARY_BUG_CONFIRMED" : "NATIVE_N_BOUNDARY_BUG_NOT_REPRODUCED");
				pp.outstream.println("NATIVE_N_BOUNDARY_DIAGNOSTIC_COMPLETE");
			}finally{DynamicDemiLog.setExponent(previous);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Invokes the inherited final method; only its add/hashAndStore destination is replaced. */
	private static Recorder measure(final Output out, final String name, final String sequence,
			final int k, final boolean quality, final int expected){
		assert(k>0 && k<32 && expected>=0) : "This fixture exercises the shared short-word loop with literal expected counts";
		final Recorder recorder=new Recorder(k, quality);
		final byte[] bases=sequence.getBytes(StandardCharsets.US_ASCII);
		final byte[] qualities=quality ? new byte[bases.length] : null;
		if(qualities!=null){Arrays.fill(qualities, (byte)40);}
		recorder.hashSmall(bases, qualities);
		write(out, name, quality ? "native_quality" : "native_no_quality", sequence, k, expected, recorder);
		return recorder;
	}

	private static void write(final Output out, final String name, final String route, final String sequence,
			final int k, final int expected, final Recorder recorder){
		assert(out!=null && recorder!=null) : "Every observed emission must reach the retained diagnostic table";
		final ByteBuilder row=new ByteBuilder();
		for(int i=0; i<Math.max(1, recorder.words.size); i++){
			row.clear().append(name).tab().append(route).tab().append(k).tab().append(sequence).tab()
				.append(expected).tab().append(recorder.words.size).tab();
			if(recorder.words.size==0){row.append("NA\tNA\tNA");}
			else{
				final long word=recorder.words.get(i);
				row.append(i).tab().append(word).tab();
				// CardinalityTracker.hashSmall uses A=0,C=1,T=2,G=3, not A/C/G/T.
				for(int base=k-1; base>=0; base--){row.append("ACTG".charAt((int)((word>>>(2*base))&3)));}
			}
			out.writer.print(row.nl());
		}
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/** Observe un-hashed canonical word emissions, with no cardinality estimate. */
	private static final class Recorder extends CardinalityTracker implements CardinalityTracker.KmerConsumer{
		Recorder(final int k, final boolean quality){super(32, k, 12345, quality ? 0.1f : 0f);}
		@Override
		public void accept(final long word, final int length){
			assert(length==k) : "The shared iterator must emit the diagnostic's declared native word length";
			add(word);
		}
		@Override
		public void hashAndStore(final long number){words.add(number);}
		@Override
		public CardinalityTracker copy(){throw new UnsupportedOperationException("Recording fixture does not copy trackers");}
		@Override
		public long cardinality(){throw new UnsupportedOperationException("Emission count is not a cardinality estimate");}
		@Override
		public void add(final CardinalityTracker other){throw new UnsupportedOperationException("Recording fixture does not merge trackers");}
		@Override
		public float[] compensationFactorLogBucketsArray(){throw new UnsupportedOperationException("No estimator calibration in emission fixture");}
		final LongList words=new LongList();
	}
}
