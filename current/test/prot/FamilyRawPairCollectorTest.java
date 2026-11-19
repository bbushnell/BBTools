package test.prot;

import ddl.DDLFamilySearch;
import prot.FamilyRawPairCollector;
import prot.FamilyRawPairCollector.Query;

/** Raw-DDL all-hit collection parity for pooled fuller negative coverage.
 * @author Yelan
 */
public final class FamilyRawPairCollectorTest{

	public static void main(String[] args){
		final String[] names={"f0","f1","f2","f3","f4"};
		final char[][] profiles={
			row(1,2,0,0),
			row(1,0,3,0),
			row(4,2,3,0),
			row(40000,0,0,0),
			row(8,8,8,0)
		};
		final char[][] weights={
			row(1,1,0,0),
			row(1,0,1,0),
			row(1,1,1,0),
			row(1,0,0,0),
			row(1,1,1,0)
		};
		final Query[] queries={
			query(0,0,1,2,3,0),
			query(0,1,1,0,0,0),
			query(1,0,4,2,3,0),
			query(2,0,40000,0,0,0),
			query(2,1,9,9,9,9)
		};
		final DDLFamilySearch.Index indexed=new DDLFamilySearch.Index(names,profiles,null,true);
		final DDLFamilySearch.Index direct=new DDLFamilySearch.Index(names,profiles,null,false);
		final DDLFamilySearch.Workspace work=new DDLFamilySearch.Workspace(names.length);
		long positives=0;
		for(int q=0; q<queries.length; q++){
			indexed.score(queries[q].sketch, work);
			positives+=work.size;
			final short[] shorts=shorts(queries[q].sketch);
			int touched=0;
			for(int f=0; f<names.length; f++){
				final int expected=direct.scoreFamily(shorts, f);
				check(work.scores[f]==expected, "Indexed raw score must match direct family scan");
				if(expected>0){touched++;}
			}
			check(work.size==touched, "Touched list must include every positive family once");
			check(work.scores[4]==0, "Untouched families must retain zero scores");
		}
		indexed.score(queries[4].sketch, work);
		check(work.size==0 && work.postings==0, "No-hit query must clear touched state");
		for(int f=0; f<names.length; f++){check(work.scores[f]==0, "No-hit query must zero stale scores");}

		final FamilyRawPairCollector.Result result=FamilyRawPairCollector.collect(indexed, queries, 2, 2);
		final FamilyRawPairCollector.FamilyStripe[] stripes=FamilyRawPairCollector.makeFamilyStripes(5, 2, 2);
		check(stripes[0].firstFamily==0 && stripes[0].lastFamily==3 &&
			stripes[1].firstFamily==3 && stripes[1].lastFamily==5,
			"Uneven family stripes must use ceil boundaries for owner() parity");
		check(stripes[FamilyRawPairCollector.owner(2, 5, 2)].owns(2),
			"Owner lookup must route a midpoint family into its stripe");
		check(FamilyRawPairCollector.owner(2, 5, 9)==2,
			"Owner lookup must clamp to the same family-owner count as stripe creation");
		final FamilyRawPairCollector.Result striped=FamilyRawPairCollector.collectStriped(indexed, queries, 2, 2, 2);
		final FamilyRawPairCollector.Result shardA=FamilyRawPairCollector.collect(indexed, queries, 0, 3, 0, queries.length, 2, 2);
		final FamilyRawPairCollector.Result shardB=FamilyRawPairCollector.collect(indexed, queries, 3, 5, 3, queries.length, 2, 2);
		final FamilyRawPairCollector.Result shallowShardA=FamilyRawPairCollector.collect(indexed, queries, 0, 3, 0, queries.length, 2, 1);
		final FamilyRawPairCollector.Result merged=FamilyRawPairCollector.merge(
			new FamilyRawPairCollector.Result[]{shardA, shardB}, names.length, queries.length, 2);
		final long[] expected={
			pair(0,0), pair(0,1), pair(0,2),
			pair(1,0), pair(1,1),
			pair(2,0), pair(2,2),
			pair(3,3)
		};
		check(result.families==names.length && result.positivePairs==positives,
			"Collector must report family count and raw positive-pair count");
		check(result.queryPairs==7 && result.familyPairs==7, "Directional top-k counts must stay separate");
		check(striped.queryPairs==result.queryPairs && striped.familyPairs==result.familyPairs &&
			striped.positivePairs==result.positivePairs && striped.size()==result.size(),
			"Family-owner striping must preserve counts and the exact serial frontier");
		final FamilyRawPairCollector.FamilyStripe[] streamStripes=FamilyRawPairCollector.makeFamilyStripes(5, 2, 2);
		final long[] streamQueries=new long[queries.length*2];
		final FamilyRawPairCollector.QueryBatch streamedBatch=FamilyRawPairCollector.collectQueryFrontier(
			indexed, queries, 0, queries.length, 0, 2, new FamilyRawPairCollector.RawPositiveObserver(){
				@Override public void observe(int family, int queryOrdinal, int hits){
					streamStripes[FamilyRawPairCollector.owner(family, names.length, streamStripes.length)].offer(family, queryOrdinal, hits);
				}
			}, streamQueries, 0);
		final int[] panels=new int[queries.length], localOrdinals=new int[queries.length];
		for(int q=0; q<queries.length; q++){
			panels[q]=queries[q].panel;
			localOrdinals[q]=queries[q].localOrdinal;
		}
		final FamilyRawPairCollector.Result streamed=FamilyRawPairCollector.finishStreaming(names.length, 2, 2,
			streamedBatch.positivePairs, streamQueries, streamedBatch.nextOffset, panels, localOrdinals, streamStripes);
		final long[] reusedQueries=new long[queries.length*2];
		final FamilyRawPairCollector.QueryBatch reused=FamilyRawPairCollector.collectQueryFrontier(
			indexed, queries, 0, queries.length, 0, 2, new DDLFamilySearch.Workspace(names.length),
			null, reusedQueries, 0);
		final long[] sliceQueries=new long[2];
		final FamilyRawPairCollector.QueryBatch slice=FamilyRawPairCollector.collectQueryFrontier(
			indexed, queries, 1, 2, 10, 2, new DDLFamilySearch.Workspace(names.length),
			null, sliceQueries, 0);
		final long[] offsetWholeQueries=new long[queries.length*2];
		final FamilyRawPairCollector.QueryBatch offsetWhole=FamilyRawPairCollector.collectQueryFrontier(
			indexed, queries, 0, queries.length, 10, 2, new DDLFamilySearch.Workspace(names.length),
			null, offsetWholeQueries, 0);
		final long[] offsetSliceQueries=new long[queries.length*2];
		final FamilyRawPairCollector.QueryBatch offsetSliceA=FamilyRawPairCollector.collectQueryFrontier(
			indexed, queries, 0, 2, 10, 2, new DDLFamilySearch.Workspace(names.length),
			null, offsetSliceQueries, 0);
		final FamilyRawPairCollector.QueryBatch offsetSliceB=FamilyRawPairCollector.collectQueryFrontier(
			indexed, queries, 2, queries.length, 12, 2, new DDLFamilySearch.Workspace(names.length),
			null, offsetSliceQueries, offsetSliceA.nextOffset);
		final int[] offsetPanels=new int[15], offsetLocalOrdinals=new int[15];
		for(int q=0; q<queries.length; q++){
			offsetPanels[10+q]=queries[q].panel;
			offsetLocalOrdinals[10+q]=queries[q].localOrdinal;
		}
		final FamilyRawPairCollector.Result offsetWholeStreamed=FamilyRawPairCollector.finishStreaming(names.length, 2, 2,
			offsetWhole.positivePairs, offsetWholeQueries, offsetWhole.nextOffset,
			offsetPanels, offsetLocalOrdinals, FamilyRawPairCollector.makeFamilyStripes(5, 2, 2));
		final FamilyRawPairCollector.Result offsetSliceStreamed=FamilyRawPairCollector.finishStreaming(names.length, 2, 2,
			offsetSliceA.positivePairs+offsetSliceB.positivePairs, offsetSliceQueries, offsetSliceB.nextOffset,
			offsetPanels, offsetLocalOrdinals, FamilyRawPairCollector.makeFamilyStripes(5, 2, 2));
		check(streamedBatch.firstOffset==0 && streamedBatch.nextOffset==result.queryPairs &&
			streamedBatch.pairs==result.queryPairs && streamedBatch.positivePairs==positives,
			"Streamed query frontier must report exact offsets and positive count");
		check(reused.nextOffset==streamedBatch.nextOffset && reused.positivePairs==streamedBatch.positivePairs,
			"Caller-owned workspace collection must match the allocating wrapper");
		for(int i=0; i<reused.nextOffset; i++){
			check(reusedQueries[i]==streamQueries[i], "Caller-owned workspace collection changed query frontier at "+i);
		}
		check(slice.nextOffset==2 && slice.positivePairs==2 &&
			sliceQueries[0]==pair(1,10) && sliceQueries[1]==pair(0,10),
			"Nonzero query slices must anchor firstQueryOrdinal to queries[from]");
		check(offsetWhole.nextOffset==offsetSliceB.nextOffset &&
			offsetWhole.positivePairs==offsetSliceA.positivePairs+offsetSliceB.positivePairs,
			"Uneven streamed slices must preserve offset pair and positive counts");
		for(int i=0; i<offsetWhole.nextOffset; i++){
			check(offsetWholeQueries[i]==offsetSliceQueries[i],
				"Uneven streamed slices changed offset query pair "+i);
		}
		check(offsetWholeStreamed.size()==offsetSliceStreamed.size(),
			"Uneven streamed slices changed shifted query-frontier union size");
		for(int i=0; i<offsetWholeStreamed.size(); i++){
			check(offsetWholeStreamed.pairAt(i)==offsetSliceStreamed.pairAt(i) &&
				offsetWholeStreamed.panelAt(i)==offsetSliceStreamed.panelAt(i) &&
				offsetWholeStreamed.localOrdinalAt(i)==offsetSliceStreamed.localOrdinalAt(i),
				"Uneven streamed slices changed shifted query-frontier panel/local pair "+i);
		}
		check(streamed.queryPairs==result.queryPairs && streamed.familyPairs==result.familyPairs &&
			streamed.positivePairs==result.positivePairs && streamed.size()==result.size(),
			"Streamed query/family union must preserve counts and the exact serial frontier");
		check(result.size()==expected.length, "Directional union must deduplicate overlapping pairs");
		for(int i=0; i<expected.length; i++){
			check(result.pairAt(i)==expected[i], "Union pair ordering changed at "+i);
			check(striped.pairAt(i)==expected[i], "Family-owner union pair ordering changed at "+i);
			check(merged.pairAt(i)==expected[i], "Shard-merged union pair ordering changed at "+i);
			check(streamed.pairAt(i)==expected[i], "Streamed union pair ordering changed at "+i);
		}
		check(merged.queryPairs==result.queryPairs && merged.familyPairs==result.familyPairs &&
			merged.positivePairs==result.positivePairs && merged.size()==result.size(),
			"Shard merge must preserve counts and the exact serial frontier");
		check(result.panelAt(7)==2 && result.localOrdinalAt(7)==0,
			"Global query ordinals must retain pooled panel/local identity");
		queries[3]=query(7,7,8,8,8,8);
		check(result.panelAt(7)==2 && result.localOrdinalAt(7)==0,
			"Result panel/local identity must not alias the caller's mutable query array");
		reject(()->FamilyRawPairCollector.collect(indexed,
			new Query[]{query(0,0,1,0,0,0),query(0,0,0,1,0,0)},1,1),
			"Duplicate panel/local queries must be rejected before pooling");
		reject(()->FamilyRawPairCollector.collectQueryFrontier(indexed,
			new Query[]{new Query(0, 0, new char[]{1,2,3})}, 0, 1, 0, 1, null, new long[1], 0),
			"Streamed raw query collection must reject wrong-width sketches explicitly");
		final FamilyRawPairCollector.FamilyStripe[] shallow=FamilyRawPairCollector.makeFamilyStripes(5, 1, 2);
		reject(()->FamilyRawPairCollector.finishStreaming(5, 2, 2, 0, new long[0], 0,
			panels, localOrdinals, shallow), "Streamed finishing must reject shallow family frontiers");
		final FamilyRawPairCollector.FamilyStripe[] partial=FamilyRawPairCollector.makeFamilyStripes(5, 2, 2);
		reject(()->FamilyRawPairCollector.finishStreaming(5, 2, 2, 0, new long[0], 0,
			panels, localOrdinals, new FamilyRawPairCollector.FamilyStripe[]{partial[0]}),
			"Streamed finishing must reject incomplete family stripe coverage");
		final FamilyRawPairCollector.FamilyStripe[] badOrdinal=FamilyRawPairCollector.makeFamilyStripes(5, 2, 2);
		badOrdinal[FamilyRawPairCollector.owner(0, 5, 2)].offer(0, 99, 1);
		reject(()->FamilyRawPairCollector.finishStreaming(5, 2, 2, 0, new long[0], 0,
			panels, localOrdinals, badOrdinal), "Streamed finishing must reject invalid family-frontier query ordinals");
		final FamilyRawPairCollector.FamilyStripe[] overdeep=FamilyRawPairCollector.makeFamilyStripes(5, 3, 2);
		overdeep[FamilyRawPairCollector.owner(0, 5, 2)].offer(0, 0, 1);
		overdeep[FamilyRawPairCollector.owner(0, 5, 2)].offer(0, 1, 1);
		overdeep[FamilyRawPairCollector.owner(0, 5, 2)].offer(0, 2, 1);
		reject(()->FamilyRawPairCollector.finishStreaming(5, 2, 2, 0, new long[0], 0,
			panels, localOrdinals, overdeep), "Streamed finishing must reject overdeep family frontiers");
		reject(()->FamilyRawPairCollector.merge(
			new FamilyRawPairCollector.Result[]{shallowShardA, shardB}, names.length, queries.length, 2),
			"Shard merge must reject family frontiers collected with insufficient local depth");
		reject(()->FamilyRawPairCollector.collect(new DDLFamilySearch.Index(names, profiles, weights, true),
			queries, 2, 2), "Weighted indexes must be rejected because raw collection ranks raw hits");
		System.out.println("FAMILY_RAW_PAIR_COLLECTOR_TEST_PASS parity=true touched=true ties=true pooled=true dedup=true");
	}

	private static Query query(int panel,int local,int a,int b,int c,int d){
		return new Query(panel, local, row(a,b,c,d));
	}

	private static char[] row(int a,int b,int c,int d){
		return new char[]{(char)a,(char)b,(char)c,(char)d};
	}

	private static short[] shorts(char[] row){
		final short[] out=new short[row.length];
		for(int i=0; i<row.length; i++){out[i]=(short)row[i];}
		return out;
	}

	private static long pair(int family,int query){
		return FamilyRawPairCollector.pack(family, query);
	}

	private static void check(boolean valid,String why){if(!valid){throw new AssertionError(why);}}

	private static void reject(Runnable action,String why){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError(why);
	}
}
