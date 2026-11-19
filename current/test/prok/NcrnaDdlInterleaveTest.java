package test.prok;

import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlInterleave;
import shared.Shared;

/** Independent literal zip/dedup orders, ties and zero-native reachability.
 * @author Brian, Raiden
 */
public final class NcrnaDdlInterleaveTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlInterleaveTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown interleave fixture parameter: "+arg);
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Pure interleave fixtures accept runtime options only");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Interleave fixtures require assertions");
			run();
			pp.outstream.println("NCRNA_DDL_INTERLEAVE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(){
		final NcrnaDdlInterleave three=new NcrnaDdlInterleave(3);
		final char[] counts={4,3,0};
		final float[] scores={3,1,2};
		three.build(new char[]{4,3,2}, scores);
		// All-positive native 0,1,2 / sketch 0,2,1 still emits 0,1,2.
		ranks(three, new int[]{1,2,3});
		three.build(counts, scores);
		// Native 0,1; sketch 0,2,1. Zip proposals 0,0,1,2,1 => 0,1,2.
		// Skipping a duplicate may not silently replace it with another proposal.
		ranks(three, new int[]{1,2,3});
		final boolean[] onlyTwo={false,false,true};
		final int[] tuple=new int[4];
		three.nativeBest(onlyTwo, tuple);
		check(Arrays.equals(tuple, new int[]{-1,-1,-1,-1}), "Native zero support must remain unreachable");
		three.sketchBest(onlyTwo, tuple);
		check(Arrays.equals(tuple, new int[]{2,2,2,2}), "The literal sketch order puts model2 second");
		three.hybridBest(onlyTwo, tuple);
		check(Arrays.equals(tuple, new int[]{3,3,3,2}), "The stable proposal zip puts model2 third without double attempts");
		check(Arrays.equals(counts, new char[]{4,3,0}) && Arrays.equals(scores, new float[]{3,1,2}),
			"Building orders borrows counts and scores read-only");

		final NcrnaDdlInterleave four=new NcrnaDdlInterleave(4);
		four.build(new char[]{8,7,0,6}, new float[]{3,2,4,1});
		// Native 0,1,3; sketch 2,0,1,3 => unique proposals 0,2,1,3.
		ranks(four, new int[]{1,3,2,4});
		final boolean[] onlyThree={false,false,false,true};
		four.nativeBest(onlyThree, tuple);
		check(Arrays.equals(tuple, new int[]{3,3,3,3}), "Positive-native rank ignores the zero-count model");
		four.hybridBest(onlyThree, tuple);
		check(Arrays.equals(tuple, new int[]{4,4,4,3}), "Unequal native/sketch lengths must still place every model once");

		four.build(new char[]{65535,65535,0,0}, new float[]{1,1,1,1});
		ranks(four, new int[]{1,2,3,4});
		final boolean[] onlyOne={false,true,false,false};
		four.nativeBest(onlyOne, tuple);
		check(Arrays.equals(tuple, new int[]{2,1,2,1}), "The maximum native count retains its complete model-index tie");
		four.sketchBest(onlyOne, tuple);
		check(Arrays.equals(tuple, new int[]{2,1,4,1}), "A sketch tie spans all four models before interleaving");
		four.hybridBest(onlyOne, tuple);
		check(Arrays.equals(tuple, new int[]{2,2,2,1}), "The hybrid is a total proposal order, not a fabricated score tie");

		four.build(new char[4], new float[]{0,-0.0f,1,0});
		ranks(four, new int[]{2,3,1,4});
		four.nativeBest(onlyOne, tuple);
		check(Arrays.equals(tuple, new int[]{-1,-1,-1,-1}), "A completely unsupported native order contributes no models");
		four.sketchBest(onlyOne, tuple);
		check(Arrays.equals(tuple, new int[]{3,2,4,1}), "Positive and negative zero must share the same model-index tie");
		four.hybridBest(onlyOne, tuple);
		check(Arrays.equals(tuple, new int[]{3,3,3,1}), "When native order is empty the hybrid becomes the complete sketch order");
		final boolean[] neither=new boolean[4];
		four.hybridBest(neither, tuple);
		check(Arrays.equals(tuple, new int[]{-1,-1,-1,-1}), "Absent oracle acceptance never invents hybrid recovery");
		final boolean[] either={true,false,true,false};
		four.hybridBest(either, tuple);
		check(Arrays.equals(tuple, new int[]{1,1,1,2}), "Multiple acceptable models select the earliest actual hybrid position");

		for(final float invalid:new float[]{-1,Float.NaN,Float.POSITIVE_INFINITY}){
			try{four.build(new char[4], new float[]{0,0,0,invalid});}
			catch(IllegalArgumentException expected){
				try{four.rank(0);}
				catch(IllegalArgumentException notReady){continue;}
				throw new AssertionError("A rejected replacement left stale hybrid ranks readable");
			}
			throw new AssertionError("Invalid sketch score defined a hybrid order: "+invalid);
		}
		largeReference();
		callerPolicy();
	}

	private static void callerPolicy(){
		final NcrnaDdlInterleave six=new NcrnaDdlInterleave(6);
		final char[] counts={4,3,1,0,0,0};
		final float[] scores={9,1,8,7,0,6};
		caller(six, counts, scores, 6, new int[]{0,1,2,3,5}, new boolean[]{false,false,true,true,true});
		caller(six, counts, scores, 3, new int[]{0,1,2}, new boolean[]{false,false,true});
		check(six.rank(3)==0 && six.rank(5)==0, "A model outside the position cap cannot retain its earlier query rank");
		final NcrnaDdlInterleave four=new NcrnaDdlInterleave(4);
		caller(four, new char[]{8,7,0,6}, new float[]{3,2,4,1}, 4,
			new int[]{0,2,1,3}, new boolean[]{false,true,false,false});
		caller(four, new char[]{0,1,1,0}, new float[]{4,0,0,3}, 4,
			new int[]{1,0,3}, new boolean[]{false,true,true});
		caller(four, new char[]{2,4,1,0}, new float[4], 4,
			new int[]{1,0}, new boolean[]{false,false});
		caller(four, new char[4], new float[]{0,2,1,-0.0f}, 4,
			new int[]{1,2}, new boolean[]{true,true});
		caller(four, new char[]{2,2,0,0}, new float[]{1,1,1,1}, 4,
			new int[]{0,1,2,3}, new boolean[]{false,false,true,true});
		caller(four, new char[4], new float[4], 4, new int[0], new boolean[0]);
		final int[] tuple=new int[4];
		four.hybridBest(new boolean[]{true,true,true,true}, tuple);
		check(Arrays.equals(tuple, new int[]{-1,-1,-1,-1}), "An empty proposal union cannot recover an acceptable but unproposed model");
		try{four.sketchFirst(0);}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("An absent model has no first-emission provenance");
	}

	private static void caller(final NcrnaDdlInterleave order, final char[] counts, final float[] scores,
			final int cap, final int[] expected, final boolean[] sketch){
		check(expected.length==sketch.length, "Each literal emitted model needs independent first-stream provenance");
		check(order.buildCaller(counts, scores, 2, cap)==expected.length && order.size()==expected.length,
			"Filtered caller policy must expose the actual union/cap length rather than a complete-permutation assumption");
		for(int position=0; position<expected.length; position++){
			final int model=expected[position];
			check(order.modelAt(position)==model && order.rank(model)==position+1 && order.sketchFirst(model)==sketch[position],
				"Literal native-first proposal zip changed its distinct position or first-emission stream at "+position);
		}
		for(int model=0; model<counts.length; model++){
			boolean present=false;
			for(final int value:expected){if(value==model){present=true;}}
			check(present || order.rank(model)==0, "An excluded model cannot retain a previous emitted position");
		}
	}

	/** An independent O(n^2) oracle selects maxima directly, without packed keys
	 * or sorting. This exercises production-size heaps, ties and zero support. */
	private static void largeReference(){
		final int models=508;
		final char[] counts=new char[models];
		final float[] scores=new float[models];
		for(int model=0; model<models; model++){
			counts[model]=(char)((model*37)%23);
			scores[model]=((model*53)%31)/7f;
		}
		final int[] expected=new int[models];
		final boolean[] usedNative=new boolean[models], usedSketch=new boolean[models], emitted=new boolean[models];
		int next=1;
		for(int round=0; round<models; round++){
			final int n=maximum(counts, usedNative);
			if(n>=0){
				usedNative[n]=true;
				if(!emitted[n]){
					emitted[n]=true;
					expected[n]=next++;
				}
			}
			final int s=maximum(scores, usedSketch);
			check(s>=0, "A complete sketch order must have a remaining model in every round");
			usedSketch[s]=true;
			if(!emitted[s]){
				emitted[s]=true;
				expected[s]=next++;
			}
		}
		check(next==models+1, "Independent repeated maximum selection must emit all original models exactly once");
		final NcrnaDdlInterleave actual=new NcrnaDdlInterleave(models);
		actual.build(counts, scores);
		ranks(actual, expected);
		for(final int cap:new int[]{50,384}){callerReference(counts, scores, cap);}
	}
	/** Independent eligible repeated maxima, then literal proposal-index zip.
	 * No packed keys, heap sorting or production emission helper is used. */
	private static void callerReference(final char[] counts, final float[] scores, final int cap){
		check(counts.length==508 && scores.length==counts.length, "The caller oracle exercises the full508-model roster");
		final int[] expected=new int[cap];
		final boolean[] sketch=new boolean[cap];
		final boolean[] usedNative=new boolean[counts.length], usedSketch=new boolean[counts.length], emitted=new boolean[counts.length];
		int eligible=0;
		for(int model=0; model<counts.length; model++){
			usedNative[model]=counts[model]<2;
			usedSketch[model]=scores[model]<=0;
			if(!usedNative[model]){eligible++;}
		}
		check(eligible>0, "This large oracle must exercise the floor; the separate literal fixture covers keep-best fallback");
		int size=0;
		for(int round=0; round<counts.length && size<cap; round++){
			final int n=maximum(counts, usedNative), s=maximum(scores, usedSketch);
			if(n<0 && s<0){break;}
			if(n>=0){
				usedNative[n]=true;
				if(!emitted[n]){
					emitted[n]=true;
					expected[size++]=n;
				}
			}
			if(s>=0 && size<cap){
				usedSketch[s]=true;
				if(!emitted[s]){
					emitted[s]=true;
					sketch[size]=true;
					expected[size++]=s;
				}
			}
		}
		check(size==cap, "The508-model oracle has enough eligible proposals to exercise the exact position cap");
		caller(new NcrnaDdlInterleave(counts.length), counts, scores, cap, expected, sketch);
	}
	private static int maximum(final char[] counts, final boolean[] used){
		assert(counts!=null && used!=null && counts.length==used.length) : "Independent native maximum selection needs one used flag per original model";
		int selected=-1;
		for(int model=0; model<counts.length; model++){
			if(!used[model] && counts[model]>0 && (selected<0 || counts[model]>counts[selected])){selected=model;}
		}
		return selected;
	}
	private static int maximum(final float[] scores, final boolean[] used){
		assert(scores!=null && used!=null && scores.length==used.length) : "Independent sketch maximum selection needs one used flag per original model";
		int selected=-1;
		for(int model=0; model<scores.length; model++){
			if(!used[model] && (selected<0 || scores[model]>scores[selected])){selected=model;}
		}
		return selected;
	}

	private static void ranks(final NcrnaDdlInterleave order, final int[] expected){
		final int[] copied=new int[expected.length];
		order.copyFirst(copied);
		for(int model=0; model<expected.length; model++){
			check(order.rank(model)==expected[model], "Literal deduplicated position differs for model "+model);
			check(copied[expected[model]-1]==model, "Copied model positions must be the independent expected permutation");
		}
		final int[] prefix=new int[Math.min(50, expected.length)];
		order.copyFirst(prefix);
		for(int position=0; position<prefix.length; position++){
			check(prefix[position]==copied[position], "Truncation must preserve positions without refilling skipped models");
		}
		copied[0]=-1;
		order.copyFirst(prefix);
		check(prefix[0]>=0, "Caller mutation of copied scratch cannot change the retained order");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
