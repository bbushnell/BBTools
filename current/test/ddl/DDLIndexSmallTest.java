package test.ddl;

import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import cardinality.DynamicDemiLog;
import ddl.DDLIndex;
import ddl.DDLIndexSmall;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/** Finite index correctness checks with literal, scalar and existing-index
 * witnesses. No external inputs, scientific ranking or timing measurements.
 * CLI handling follows the native template/finite-fixture lifecycle.
 * @author Brian Bushnell, Raiden
 */
public final class DDLIndexSmallTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DDLIndexSmallTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown index-test parameter: "+arg);
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& parser.maxReads<0 && parser.samplerate==1, "Index fixtures use complete synthetic arrays without files or sampling");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Index boundary and overflow witnesses require assertions enabled");
			literal();
			nativeObjects();
			widths();
			ordinals();
			invalid();
			pp.outstream.println("DDL_INDEX_SMALL_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static void literal(){
		final char[][] references={{65535, 0, 7}, {65535, 32768, 0}, {32768, 0, 0}, {0, 1, 7}};
		final char[] query={65535, 32768, 7};
		final DDLIndexSmall small=new DDLIndexSmall(references);
		final DDLIndex matrix=new DDLIndex(3);
		for(int id=0; id<references.length; id++){
			int filled=0;
			for(final char key:references[id]){if(key!=0){filled++;}}
			matrix.add(id, references[id], filled);
		}
		final int[] expected={2, 2, 0, 1};
		check(Arrays.equals(small.query(query), expected) && Arrays.equals(matrix.query(query), expected),
			"Literal signed-boundary keys must agree with the existing matrix DDL index");
		check(small.numClades()==4 && small.buckets()==3 && small.populatedBuckets()==3
			&& small.emptyBuckets()==0 && small.pointerSlots()==3 && small.postingPayloadBytes()==28
			&& small.populatedCells()==5 && matrix.populatedCells()==5,
			"Seven postings occupy five distinct bin/key cells in three populated bins");
		final int[] counts={9, 9, 9, 9};
		check(small.query(query, counts)==counts && Arrays.equals(counts, expected),
			"Query reuses and replaces the caller's complete counter array");
		check(small.accumulateCounts(query, counts)==5 && Arrays.equals(counts, new int[]{4, 4, 0, 2}),
			"Additive lookup retains prior counts and reports exactly the five visited postings");
		check(small.accumulateCounts(new char[3], counts)==0 && counts[0]==4,
			"Empty bins do not match or clear existing additive counts");
		small.query(new char[3], counts);
		check(Arrays.equals(counts, new int[4]), "An empty query replacement clears all previous counts");
		final int[] raw={65535, 32768, (1<<16)|1, (1<<16)|32768, (2<<16)|7, 0};
		check(Arrays.equals(small.query(raw, 5), new int[]{2, 2, 1, 2}),
			"Distinct raw keys in one bin may match different references, while unused capacity is ignored");
		small.query(raw, 3, counts);
		check(Arrays.equals(counts, new int[]{1, 1, 1, 1}), "Only the declared live query prefix contributes");
		final int[][] parts={new int[4], new int[4], new int[4]};
		check(small.query(raw, 5, new int[]{1, 2, 3}, parts)==parts
			&& Arrays.equals(parts[0], new int[]{1, 1, 1, 0})
			&& Arrays.equals(parts[1], new int[]{0, 1, 0, 1})
			&& Arrays.equals(parts[2], new int[]{1, 0, 0, 1}),
			"Exclusive bin ends assign each raw match to its declared partition");
		check(small.accumulateCounts(raw, 5, new int[]{1, 2, 3}, parts)==7
			&& Arrays.equals(parts[2], new int[]{2, 0, 0, 2}),
			"Partition accumulation adds instead of replacing and reports actual postings visited");
		small.query(raw, 0, new int[]{1, 2, 3}, parts);
		for(final int[] row:parts){check(Arrays.equals(row, new int[4]), "Empty partition queries replace previous counters");}
		references[0][0]=0;
		references[1]=new char[3];
		check(Arrays.equals(small.query(query), expected), "Producer array edits cannot alter the index's copied postings");
		final DDLIndexSmall empty=new DDLIndexSmall(new char[4][3]);
		check(empty.populatedBuckets()==0 && empty.populatedCells()==0 && empty.emptyBuckets()==3
			&& empty.postingPayloadBytes()==0 && Arrays.equals(empty.query(query), new int[4]),
			"All-empty references retain their ordinals without manufacturing empty/empty hits");
	}

	/** Independent O(query keys times references) oracle, including arbitrary
	 * and maximum bin counts. No boxed keys or per-query-key objects are created. */
	private static void widths(){
		final Random random=new Random(20261114L);
		for(final int width:new int[]{1, 3, 85, 170, 256, 65536}){
			final int models=width==65536 ? 7 : 23;
			final char[][] references=new char[models][width];
			for(final char[] row:references){
				for(int bin=0; bin<width; bin++){
					if(random.nextInt(5)==0){row[bin]=(char)(1+random.nextInt(65535));}
				}
			}
			final char[] query=new char[width];
			final int[] raw=new int[2*width];
			for(int bin=0; bin<width; bin++){
				final int key=bin%3==0 && references[0][bin]!=0 ? references[0][bin] : 1+random.nextInt(65535);
				query[bin]=bin%7==0 ? 0 : (char)key;
				raw[2*bin]=(bin<<16)|key;
				raw[2*bin+1]=(bin<<16)|(1+key%65535);
			}
			final int[] expected=new int[models], expectedRaw=new int[models];
			for(int id=0; id<models; id++){
				for(int bin=0; bin<width; bin++){
					if(query[bin]!=0 && query[bin]==references[id][bin]){expected[id]++;}
				}
			}
			for(final int packed:raw){
				for(int id=0; id<models; id++){
					if(references[id][packed>>>16]==(packed&65535)){expectedRaw[id]++;}
				}
			}
			final DDLIndexSmall index=new DDLIndexSmall(references);
			check(Arrays.equals(index.query(query), expected) && Arrays.equals(index.query(raw, raw.length), expectedRaw),
				"Indexed counts must agree with independent scalar membership at width="+width);
			final int[] scratch=new int[models];
			for(int repeat=0; repeat<3; repeat++){
				index.query(raw, raw.length, scratch);
				check(Arrays.equals(scratch, expectedRaw), "Repeated raw replacements must not accumulate old counts");
				index.query(query, scratch);
				check(Arrays.equals(scratch, expected), "Raw and query-DDL calls may safely share worker scratch");
			}
		}
	}

	private static void ordinals(){
		final char[][] references=new char[65536][1];
		references[0][0]=1;
		references[32768][0]=32768;
		references[65535][0]=65535;
		final DDLIndexSmall index=new DDLIndexSmall(references);
		final int[] counts=index.query(new int[]{1, 32768, 65535}, 3);
		check(counts[0]==1 && counts[32768]==1 && counts[65535]==1 && index.postingPayloadBytes()==12,
			"Unsigned reference65535 and key65535 must survive the packed posting value minus one");
		int total=0;
		for(final int count:counts){total+=count;}
		check(total==3, "Maximum ordinal lookup cannot introduce matches in other reference slots");
	}

	private static void invalid(){
		reject(()->new DDLIndexSmall((char[][])null));
		reject(()->new DDLIndexSmall(new char[0][]));
		reject(()->new DDLIndexSmall(new char[][]{null}));
		reject(()->new DDLIndexSmall(new char[][]{{1}, {1, 2}}));
		reject(()->new DDLIndexSmall(new char[][]{new char[0]}));
		reject(()->new DDLIndexSmall(new char[][]{new char[65537]}));
		final DDLIndexSmall index=new DDLIndexSmall(new char[][]{{1, 2}, {3, 4}});
		reject(()->index.query(new char[1]));
		reject(()->index.query(new char[2], new int[1]));
		reject(()->index.query(new int[]{1}, -1));
		reject(()->index.query(new int[]{1}, 2));
		reject(()->index.query(new int[]{0}, 1));
		reject(()->index.query(new int[]{(2<<16)|1}, 1));
		final int[] alias={1, 2};
		reject(()->index.query(alias, 2, alias));
		reject(()->index.query(new int[]{1}, 1, new int[]{1}, new int[][]{new int[2]}));
		reject(()->index.query(new int[]{1}, 1, new int[]{1, 1, 2}, new int[3][2]));
		reject(()->index.query(new int[]{1}, 1, new int[]{1, 2}, new int[][]{alias, alias}));
		boolean overflow=false;
		try{index.accumulateCounts(new char[]{1, 0}, new int[]{Integer.MAX_VALUE, 0});}
		catch(AssertionError expected){overflow=true;}
		check(overflow, "Additive scratch overflow must fail loudly rather than publish wrapped match counts");
	}

	/** Native objects retain original ordinals, unsigned keys, empty-bin rules
	 * and caller scratch; index construction owns its postings after reuse. */
	private static void nativeObjects(){
		final int previous=DynamicDemiLog.exponentBits();
		try{
			DynamicDemiLog.setExponent(4);
			final char[][] keys={{65535,0,7,0}, {65535,32768,0,0}, {32768,0,0,0}, {0,1,7,0}};
			final DynamicDemiLog[] references=new DynamicDemiLog[keys.length];
			for(int id=0; id<keys.length; id++){references[id]=DynamicDemiLog.fromArray(keys[id], id, null, 7, -1);}
			final DynamicDemiLog query=DynamicDemiLog.fromArray(new char[]{65535,32768,7,0}, -1, null, 7, -1);
			final DDLIndexSmall index=new DDLIndexSmall(references);
			final int[] expected={2,2,0,1}, scratch={9,9,9,9}, direct=new int[4];
			check(Arrays.equals(index.query(query), expected), "Native query results must retain all original model ordinals, including zeros");
			check(index.query(query, scratch)==scratch && Arrays.equals(scratch, expected),
				"Native query scratch is replaced without a result allocation");
			for(int id=0; id<references.length; id++){
				DynamicDemiLog.compareDetailed(query.maxArray(), references[id].maxArray(), direct);
				check(scratch[id]==direct[1], "Every native object score must equal the shared nonempty equal-bin count");
			}
			check(index.accumulateCounts(query, scratch)==5 && Arrays.equals(scratch, new int[]{4,4,0,2}),
				"The native additive overload must retain existing counters and report actual postings visited");
			references[0].clear();
			references[1]=null;
			check(Arrays.equals(index.query(query), expected), "Source object reuse/replacement cannot mutate copied postings");
			query.clear();
			index.query(query, scratch);
			check(Arrays.equals(scratch, new int[4]), "Empty native queries clear stale counters without empty-empty matches");
			reject(()->new DDLIndexSmall((DynamicDemiLog[])null));
			reject(()->new DDLIndexSmall(new DynamicDemiLog[0]));
			reject(()->new DDLIndexSmall(new DynamicDemiLog[]{null}));
			reject(()->new DDLIndexSmall(new DynamicDemiLog[]{query, DynamicDemiLog.create(8, 7, 12345, 0, false)}));
			reject(()->index.query((DynamicDemiLog)null));
			reject(()->index.query(DynamicDemiLog.create(8, 7, 12345, 0, false), scratch));
			reject(()->index.query(query, new int[3]));
			final DynamicDemiLog actual=DynamicDemiLog.create(128, 7, 12345, 0, false);
			actual.hash("ACGTAGCTAGGCTAGCTAGCATGCATGCTAGCTAGCTAGCATG".getBytes(java.nio.charset.StandardCharsets.US_ASCII), null);
			final DDLIndexSmall actualIndex=new DDLIndexSmall(new DynamicDemiLog[]{actual});
			check(actualIndex.query(actual)[0]==actual.filledBuckets(), "Native factory/insertion objects must match all of their nonempty bins");
		}finally{DynamicDemiLog.setExponent(previous);}
	}

	private static void reject(final Runnable operation){
		boolean rejected=false;
		try{operation.run();}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Invalid dimensions, live keys or scratch aliases must fail for their explicit contract");
	}

	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
