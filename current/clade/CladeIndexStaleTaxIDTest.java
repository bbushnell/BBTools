package clade;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;

import bin.AdjustEntropy;
import cardinality.DynamicDemiLog;
import ddl.DDLIndexBase;
import ddl.DDLRecord;
import structures.IntHashMap;
import tax.TaxNode;
import tax.TaxTree;

/** Regression fixtures for sketch-LCA fallbacks when local DDL records outpace taxtree.tsv.gz.
 * @author Brian Bushnell, Yelan */
public final class CladeIndexStaleTaxIDTest {

	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {AdjustEntropy.load(4, 150);}
		final TaxTree tree=tree();
		final TaxTree oldShared=TaxTree.getTree(), oldClade=CladeObject.tree;
		try{
			setSharedTree(tree);
			CladeObject.tree=tree;
			check(CladeIndex.treeContains(tree, 2), "fixture tree should contain domain taxid 2");
			check(!CladeIndex.treeContains(tree, 5), "fixture tree should miss in-range deleted taxid");
			check(!CladeIndex.treeContains(tree, 1341699), "fixture tree should miss stale sketch taxid");
			Clade staleA=clade(1341699, TaxTree.SPECIES, "d__Bacteria;p__Actinomycetota;c__Actinomycetes");
			Clade staleB=clade(42, -1, "d__Bacteria;p__Actinomycetota;c__Acidimicrobiia");
			check(CladeIndex.sketchLCA(staleA, staleA.taxID, staleA.level, staleA, tree)==TaxTree.SPECIES,
				"exact stale taxids should use stored clade level");
			check(CladeIndex.sketchLCA(staleA, staleB.taxID, staleA.level, staleB, tree)==TaxTree.PHYLUM,
				"stale local sketch taxids should fall back to lineage text instead of tree LCA");
			check(CladeIndex.sketchLCA(clade(2, TaxTree.SPECIES, null), 2, TaxTree.SPECIES,
				clade(2, TaxTree.SPECIES, null), tree)==TaxTree.DOMAIN,
				"valid taxids should use the active tree rank before stored Clade.level");
			check(CladeIndex.sketchLCA(clade(4, TaxTree.SPECIES, null), 3, TaxTree.SPECIES,
				clade(3, TaxTree.SPECIES, null), tree)==TaxTree.DOMAIN,
				"merged taxids should resolve through the active tree before stored Clade.level");
			check(CladeIndex.sketchLCA(clade(2, -1, null), 3, TaxTree.SPECIES,
				clade(3, -1, null), tree)==TaxTree.LIFE,
				"valid taxids should use the taxonomy tree LCA");
			check(CladeIndex.sketchLCA(clade(5, TaxTree.SPECIES, null), 5, -1,
				clade(5, -1, null), tree)==-1, "in-range deleted taxids should be unknown-safe");
			check(CladeIndex.sketchLCA(clade(1341699, TaxTree.SPECIES, null), 42, -1,
				clade(42, -1, null), tree)==-1, "out-of-range missing cached lineages should be unknown-safe");
			Comparison comparison=new Comparison();
			comparison.ref=clade(1341699, -1, null);
			check(QueryResult.lcaFor(comparison, null, clade(42, -1, null))==TaxTree.LIFE,
				"QueryResult lineage fallback should not assert on missing cached lineages");
			checkAddSketchInfoPath(tree);
			checkAddSketchInfoSharedTreeIsolation();
			checkQueryResultPrivateSketchCap();
			checkDownstreamConsumers();
			System.out.println("CladeIndexStaleTaxIDTest PASS");
		}finally{
			setSharedTree(oldShared);
			CladeObject.tree=oldClade;
		}
	}

	private static void checkAddSketchInfoPath(TaxTree tree){
		final DynamicDemiLog ddl=DynamicDemiLog.create(4, 5, 1L, 0f, true);
		ddl.hashAndStore(123);
		final ArrayList<DDLRecord> records=new ArrayList<DDLRecord>();
		records.add(new DDLRecord(ddl, 1, 1341699, "stale"));
		records.add(new DDLRecord(ddl, 2, 42, "staleOther"));
		final CladeIndex index=new CladeIndex(new ArrayList<Clade>());
		index.ddlIndex=DDLIndexBase.create(4);
		index.ddlIndex.addAll(records, 1);
		index.sketchRecords=records;
		index.cladeMap=new ConcurrentHashMap<Integer, Clade>();
		index.cladeMap.put(42, clade(42, -1, "d__Bacteria;p__Actinomycetota;c__Acidimicrobiia"));
		final Clade query=clade(0, -1, null);
		query.ddl=ddl;
		final Comparison comp=new Comparison();
		comp.query=query;
		comp.ref=clade(1341699, TaxTree.SPECIES, null);
		final ArrayList<Comparison> results=new ArrayList<Comparison>();
		results.add(comp);
		final int oldMin=CladeIndex.minSketchMatches, oldMax=CladeIndex.maxSketchHits;
		try{
			CladeIndex.minSketchMatches=0;
			CladeIndex.maxSketchHits=2;
			index.addSketchInfo(results, query);
		}finally{
			CladeIndex.minSketchMatches=oldMin;
			CladeIndex.maxSketchHits=oldMax;
		}
		check(comp.sketchLCA==TaxTree.SPECIES, "addSketchInfo should survive exact stale taxids");
		check(results.size()==2 && results.get(1).isSketchHit && results.get(1).sketchLCA<0,
			"addSketchInfo should safely add a sketch-only stale hit");
	}

	private static void checkAddSketchInfoSharedTreeIsolation(){
		final DynamicDemiLog ddl=DynamicDemiLog.create(4, 5, 1L, 0f, true);
		ddl.hashAndStore(123);
		final ArrayList<DDLRecord> records=new ArrayList<DDLRecord>();
		records.add(new DDLRecord(ddl, 1, 2, "loadedTreeTaxid"));
		final CladeIndex index=new CladeIndex(new ArrayList<Clade>());
		index.ddlIndex=DDLIndexBase.create(4);
		index.ddlIndex.addAll(records, 1);
		index.sketchRecords=records;
		index.cladeMap=new ConcurrentHashMap<Integer, Clade>();
		index.cladeMap.put(2, clade(2, TaxTree.DOMAIN, "d__Bacteria;p__Actinomycetota"));
		final Clade query=clade(0, -1, null);
		query.ddl=ddl;
		final Comparison comp=new Comparison();
		comp.query=query;
		comp.ref=clade(3, TaxTree.DOMAIN, "d__Bacteria;p__Actinomycetota");
		final ArrayList<Comparison> results=new ArrayList<Comparison>();
		results.add(comp);
		final int oldMin=CladeIndex.minSketchMatches;
		try{
			CladeIndex.minSketchMatches=0;
			index.addSketchInfo(results, query, 1, false);
		}finally{
			CladeIndex.minSketchMatches=oldMin;
		}
		check(comp.sketchLCA==TaxTree.PHYLUM,
			"addSketchInfo false should ignore shared-tree taxid LCA and use cached lineage text");
	}

	private static void checkQueryResultPrivateSketchCap(){
		final Clade query=clade(0, -1, null);
		final ArrayList<Comparison> hits=new ArrayList<Comparison>();
		for(int i=0; i<10; i++){
			Comparison c=comparison(query, clade(100+i, TaxTree.SPECIES, null), 0);
			c.isSketchHit=true;
			c.wkid=1f;
			c.kmerMatches=10-i;
			hits.add(c);
		}
		final int oldMax=CladeIndex.maxSketchHits;
		final boolean oldDDL=Clade.MAKE_DDLS;
		try{
			CladeIndex.maxSketchHits=5;
			Clade.MAKE_DDLS=false;
			final QueryResult result=QueryResult.build(query, new FixedSketchIndex(hits),
				1, 10, 1, false, null, 10);
			check(result.displayList.size()==10,
				"QueryResult should retain the per-call sketch cap, not global maxSketchHits");
		}finally{
			CladeIndex.maxSketchHits=oldMax;
			Clade.MAKE_DDLS=oldDDL;
		}
	}

	private static void checkDownstreamConsumers(){
		if(!CladeRanking.ready()){throw new AssertionError("ranking.bbnet was not loaded");}
		if(!CladeConfidence.v2Ready()){throw new AssertionError("confidence.bbnets.gz was not loaded");}
		final Clade query=clade(0, -1, null);
		query.setBases(1000000);
		final ArrayList<Comparison> group=new ArrayList<Comparison>();
		group.add(comparison(query, clade(1341699, -1, null), 2));
		group.add(comparison(query, clade(5, -1, null), 1));
		CladeRanking.scoreAndResort(group, 32768);
		for(int i=0; i<group.size(); i++){
			group.get(i).cacheConfidence(group, i);
			group.get(i).toBytes(null);
		}
	}

	private static Comparison comparison(Clade query, Clade ref, float composite){
		Comparison c=new Comparison();
		c.query=query;
		c.ref=ref;
		c.composite=composite;
		c.gcdif=c.strdif=c.hhdif=c.cagadif=c.k3dif=c.k4dif=c.k5dif=0.01f;
		c.ssudif=1;
		c.kid=c.wkid=0.99f;
		c.kmerMatches=100;
		c.sketchLCA=-1;
		return c;
	}

	private static Clade clade(int tid, int level, String lineage){
		Clade c=new Clade(tid, level, "c"+tid);
		c.lineage=lineage;
		c.add("ACGTACGTACGTACGTACGT".getBytes(), null);
		c.finish();
		return c;
	}

	private static TaxTree tree() throws Exception{
		final TaxNode[] nodes=new TaxNode[6];
		nodes[1]=new TaxNode(1, 1, TaxTree.LIFE, TaxTree.LIFE, "root");
		nodes[2]=new TaxNode(2, 1, TaxTree.DOMAIN, TaxTree.DOMAIN, "Bacteria");
		nodes[3]=new TaxNode(3, 1, TaxTree.DOMAIN, TaxTree.DOMAIN, "Archaea");
		final IntHashMap merged=new IntHashMap();
		merged.put(4, 3);
		final Constructor<TaxTree> c=TaxTree.class.getDeclaredConstructor(TaxNode[].class,
			IntHashMap.class, int.class, boolean.class, boolean.class, boolean.class, int.class);
		c.setAccessible(true);
		return c.newInstance(nodes, merged, 1, false, false, false, TaxTree.DOMAIN);
	}

	private static void setSharedTree(TaxTree tree) throws Exception{
		Field field=TaxTree.class.getDeclaredField("sharedTree");
		field.setAccessible(true);
		field.set(null, tree);
	}

	private static void check(boolean ok, String msg){
		if(!ok){throw new AssertionError(msg);}
	}

	private static final class FixedSketchIndex extends CladeIndex {

		FixedSketchIndex(ArrayList<Comparison> hits_){
			super(new ArrayList<Clade>());
			hits=hits_;
		}

		@Override
		public ArrayList<Comparison> findBest(final Clade c, final int maxHits, final int maxSketchHits){
			return new ArrayList<Comparison>(hits);
		}

		private final ArrayList<Comparison> hits;
	}
}
