package clade;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import parse.Parse;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Audit tool: compares the TREE LCA (tree.commonAncestorLevel) against the STRING LCA
 * (CladeIndex.lineageLCA) on real reference (taxID, lineage) pairs from a spectra DB.
 * These are the two paths CladeIndex.addSketchInfo takes for sketchLCA in usetree=t vs
 * usetree=f mode. The sketch anchor is tree-INDEPENDENT, so the per-pair divergence here is
 * the raw sketchLCA-level divergence between the modes.
 *
 * ★ HOW THAT DIVERGENCE PROPAGATES (do NOT read encodeLevel9 agreement as "what the net sees"):
 * The RAW sketchLCA level feeds Comparison.compositeScore via lcaFactor=(18-lca)^2 (289 at
 * SUBSPECIES(1) vs 256 at SPECIES(2)), and `composite` then feeds BOTH the pre-ranking sort AND
 * the ranking net's `norm`/`spread` dims (CladeRanking.buildVector). CladeRanking's OWN sketchLCA
 * dim, encodeLevel9(level), CLAMPS SUBSPECIES(1)==SPECIES(2)==1.0 -- so encodeLevel9 agreement
 * MASKS the very divergence that flips top hits. Track BOTH: raw-level parity (drives compositeScore)
 * is the load-bearing one; encodeLevel9 parity is a weaker, misleading proxy.
 *
 * Reports: level-exact agreement, encodeLevel9 agreement (labeled as the WEAK proxy it is), a
 * (strLevel x treeLevel) confusion matrix, and labeled sample lineage pairs per divergence class
 * -- so the CAUSE of each mismatch is visible, not inferred.
 *
 * @author Noire
 */
public class LcaAudit {

	static boolean selfAll=false;//iterate ALL records: lineageLCA(self,self) vs commonAncestorLevel(tid,tid)
	static boolean useFresh=false;//sample pairs on FRESH lineages (Clade.lineage, with the fix) not stored DB ones

	public static void main(String[] args){
		String spectra=null, treePath=null, nodes=null;
		long pairs=2000000;
		long seed=17;
		int maxSamplesPerClass=8;
		for(String a : args){
			String[] s=a.split("=", 2);
			String k=s[0].toLowerCase();
			String v=s.length>1 ? s[1] : null;
			if(k.equals("spectra") || k.equals("in") || k.equals("ref")){spectra=v;}
			else if(k.equals("tree")){treePath=v;}
			else if(k.equals("pairs")){pairs=Parse.parseKMG(v);}
			else if(k.equals("seed")){seed=Long.parseLong(v);}
			else if(k.equals("samples")){maxSamplesPerClass=Integer.parseInt(v);}
			else if(k.equals("nodes")){nodes=v;}
			else if(k.equals("selfall")){selfAll=Parse.parseBoolean(v);}
			else if(k.equals("fresh")){useFresh=Parse.parseBoolean(v);}
			else{assert(false) : "Unknown arg "+a;}
		}
		assert(spectra!=null && treePath!=null) : "Need spectra= and tree=";

		System.err.println("Loading tree from "+treePath);
		TaxTree tree=TaxTree.loadTaxTree(treePath, System.err, false, false);
		CladeObject.tree=tree;//so Clade.lineage(tid) can regenerate lineages fresh from the tree (with the fix)

		System.err.println("Loading clades from "+spectra);
		Clade.MAKE_DDLS=false;
		CladeLoader loader=new CladeLoader();
		ConcurrentHashMap<Integer, Clade> map=loader.load(java.util.Arrays.asList(spectra), null);

		ArrayList<Clade> valid=new ArrayList<Clade>();
		int noLineage=0, notInTree=0;
		for(Clade c : map.values()){
			if(c.taxID<1){continue;}
			if(c.lineage==null){noLineage++; continue;}
			if(tree.getNode(c.taxID, true)==null){notInTree++; continue;}
			valid.add(c);
		}
		System.err.println("Loaded "+map.size()+" clades; valid="+valid.size()
			+", noLineage="+noLineage+", notInTree="+notInTree);
		if(valid.size()<2){System.err.println("Not enough valid clades."); return;}

		if(nodes!=null){dumpNodes(tree, map, nodes); return;}

		if(selfAll){selfCheckAll(tree, valid); return;}

		//Bucket by family token for NEAR sampling (fine-level LCAs, where ranking is sensitive).
		HashMap<String, ArrayList<Clade>> famMap=new HashMap<String, ArrayList<Clade>>();
		for(Clade c : valid){
			String fam=token(c.lineage, ";f__");
			if(fam==null){continue;}
			ArrayList<Clade> list=famMap.get(fam);
			if(list==null){famMap.put(fam, list=new ArrayList<Clade>());}
			list.add(c);
		}
		ArrayList<ArrayList<Clade>> famLists=new ArrayList<ArrayList<Clade>>();
		for(ArrayList<Clade> list : famMap.values()){if(list.size()>=2){famLists.add(list);}}
		System.err.println("Families with >=2 members: "+famLists.size());

		Random rng=new Random(seed);
		final int n=valid.size();

		//Tally: confusion[strIdx][treeIdx], where index = level+1 (so -1 maps to 0). Levels 0..11 -> 1..12.
		long[][] confusion=new long[13][13];
		long total=0, agreeLevel=0, agreeEnc=0;
		double sumAbsEncDiff=0;
		//Per divergence class (strLevel,treeLevel), keep samples.
		HashMap<String, ArrayList<String>> samples=new HashMap<String, ArrayList<String>>();
		HashMap<String, Long> classCount=new HashMap<String, Long>();

		long half=pairs/2;
		for(long p=0; p<pairs; p++){
			Clade a, b;
			if(p<half && !famLists.isEmpty()){//NEAR pair: two distinct members of one family
				ArrayList<Clade> fam=famLists.get(rng.nextInt(famLists.size()));
				a=fam.get(rng.nextInt(fam.size()));
				b=fam.get(rng.nextInt(fam.size()));
				int guard=0;
				while(b.taxID==a.taxID && guard++<8){b=fam.get(rng.nextInt(fam.size()));}
				if(b.taxID==a.taxID){continue;}
			}else{//FAR pair: two random valid clades
				a=valid.get(rng.nextInt(n));
				b=valid.get(rng.nextInt(n));
				if(a.taxID==b.taxID){continue;}
			}

			int treeLevel=tree.commonAncestorLevel(a.taxID, b.taxID);
			CharSequence linA=useFresh ? Clade.lineage(a.taxID) : a.lineage;
			CharSequence linB=useFresh ? Clade.lineage(b.taxID) : b.lineage;
			int strLevel=CladeIndex.lineageLCA(linA, linB);
			float encTree=enc(treeLevel);
			float encStr=enc(strLevel);

			total++;
			confusion[idx(strLevel)][idx(treeLevel)]++;
			if(treeLevel==strLevel){agreeLevel++;}
			float d=Math.abs(encTree-encStr);
			if(d<1e-6f){agreeEnc++;}
			else{
				sumAbsEncDiff+=d;
				String key=strLevel+"->"+treeLevel;//string returned strLevel, tree returned treeLevel
				classCount.merge(key, 1L, Long::sum);
				ArrayList<String> list=samples.get(key);
				if(list==null){samples.put(key, list=new ArrayList<String>());}
				if(list.size()<maxSamplesPerClass){
					list.add("  tid "+a.taxID+" ["+lname(tree,a.taxID)+"]  vs  tid "+b.taxID+" ["+lname(tree,b.taxID)+"]"
						+"\n    A: "+a.lineage
						+"\n    B: "+b.lineage);
				}
			}
		}

		System.out.println("===== LCA AUDIT (tree.commonAncestorLevel vs CladeIndex.lineageLCA) =====");
		System.out.println("pairs sampled: "+total+"  (seed="+seed+")");
		System.out.println("level-exact agreement:   "+agreeLevel+"/"+total+" = "+fmt(agreeLevel/(double)total));
		System.out.println("encodeLevel9 agreement:  "+agreeEnc+"/"+total+" = "+fmt(agreeEnc/(double)total)
			+"   (this is what the ranking net actually sees)");
		System.out.println("mean |encTree-encStr| over ALL pairs: "+fmt(sumAbsEncDiff/total));
		System.out.println();

		//Confusion matrix (only nonempty rows/cols), levels labeled.
		System.out.println("Confusion matrix rows=strLevel(lineageLCA), cols=treeLevel(commonAncestorLevel):");
		System.out.print("str\\tree");
		for(int t=0; t<13; t++){if(colUsed(confusion,t)){System.out.print("\t"+lvlName(t-1));}}
		System.out.println();
		for(int r=0; r<13; r++){
			boolean rowUsed=false;
			for(int t=0; t<13; t++){if(confusion[r][t]>0){rowUsed=true; break;}}
			if(!rowUsed){continue;}
			System.out.print(lvlName(r-1));
			for(int t=0; t<13; t++){if(colUsed(confusion,t)){System.out.print("\t"+confusion[r][t]);}}
			System.out.println();
		}
		System.out.println();

		//Divergence classes sorted by count desc.
		System.out.println("===== ENCODE-DIVERGENT CLASSES (strLevel -> treeLevel), by frequency =====");
		ArrayList<String> keys=new ArrayList<String>(classCount.keySet());
		keys.sort((x,y)->Long.compare(classCount.get(y), classCount.get(x)));
		for(String key : keys){
			String[] st=key.split("->");
			int sl=Integer.parseInt(st[0]), tl=Integer.parseInt(st[1]);
			System.out.println("\n["+lvlName(sl)+" (str) -> "+lvlName(tl)+" (tree)]  count="+classCount.get(key)
				+"   encStr="+fmt(enc(sl))+" encTree="+fmt(enc(tl))+" |diff|="+fmt(Math.abs(enc(sl)-enc(tl))));
			for(String samp : samples.get(key)){System.out.println(samp);}
		}
	}

	/** Dump the ancestor chain + level semantics for specific taxIDs (comma-separated), plus pairwise LCAs. */
	private static void dumpNodes(TaxTree tree, ConcurrentHashMap<Integer, Clade> map, String csv){
		String[] ids=csv.split(",");
		for(String idStr : ids){
			int id=Integer.parseInt(idStr.trim());
			TaxNode tn=tree.getNode(id, true);
			System.out.println("\n==== taxID "+id+" ====");
			if(tn==null){System.out.println("  NOT IN TREE"); continue;}
			Clade c=map.get(id);
			System.out.println("  stored lineage (DB, pre-fix): "+(c!=null ? c.lineage : "<no clade in DB>"));
			System.out.println("  FRESH lineage (Clade.lineage, with fix): "+Clade.lineage(id));
			System.out.println("  ancestor chain (leaf -> species):");
			TaxNode t=tn; int guard=0;
			while(t!=null && guard++<20){
				System.out.println(String.format("    id=%-9d level=%-12s levelE=%-16s isSimple=%-5s name=%s",
					t.id, lvlName(t.level), TaxTree.levelToStringExtended(t.levelExtended),
					TaxTree.isSimple(t.levelExtended), t.name));
				if(t.level>=TaxTree.SPECIES || t.id==t.pid){break;}
				t=tree.getNode(t.pid);
			}
		}
		//SELF-LCA is the case that drives the flip: sketchLCA = LCA(ref, best-sketch-anchor), and for a strain
		//reference the anchor is that same strain. Uses FRESH lineages (Clade.lineage, with the fix) so this
		//verifies the fix WITHOUT regenerating the DB. tree.commonAncestorLevel(x,x)=x's formal level.
		System.out.println("\n==== SELF-LCA (fresh lineage, THE flip case) ====");
		for(String idStr : ids){
			int a=Integer.parseInt(idStr.trim());
			CharSequence lin=Clade.lineage(a);
			int treeL=tree.commonAncestorLevel(a, a);
			int strL=CladeIndex.lineageLCA(lin, lin);
			String flag=(treeL==strL ? "MATCH" : "*** MISMATCH ***");
			System.out.println("  "+a+" self: commonAncestorLevel="+lvlName(treeL)+"  lineageLCA="+lvlName(strL)+"  "+flag);
		}
		System.out.println("\n==== pairwise (commonAncestorLevel vs lineageLCA, FRESH lineages) ====");
		for(int i=0; i<ids.length; i++){
			for(int j=i+1; j<ids.length; j++){
				int a=Integer.parseInt(ids[i].trim()), b=Integer.parseInt(ids[j].trim());
				int treeL=tree.commonAncestorLevel(a, b);
				int strL=CladeIndex.lineageLCA(Clade.lineage(a), Clade.lineage(b));
				String flag=(treeL==strL ? "MATCH" : "*** MISMATCH ***");
				System.out.println("  "+a+" vs "+b+": commonAncestorLevel="+lvlName(treeL)+"  lineageLCA="+lvlName(strL)+"  "+flag);
			}
		}
	}

	/** Whole-DB self-LCA parity: for EVERY record, does lineageLCA(freshLineage,freshLineage) equal
	 *  tree.commonAncestorLevel(tid,tid)? This is the sketchLCA self case (ref==anchor) generalized to all
	 *  records -- the comprehensive validation of the emitter+lineageLCA fix, independent of the 8000-shred A/B. */
	private static void selfCheckAll(TaxTree tree, ArrayList<Clade> valid){
		long total=0, match=0, fastMatch=0;
		HashMap<String, Long> mism=new HashMap<String, Long>();
		HashMap<String, ArrayList<String>> samp=new HashMap<String, ArrayList<String>>();
		for(Clade c : valid){
			CharSequence lin=Clade.lineage(c.taxID);
			int treeL=tree.commonAncestorLevel(c.taxID, c.taxID);
			int strL=CladeIndex.lineageLCA(lin, lin);
			//Fast-path (runtime) self-LCA: exact-taxID branch uses stored level when >=0, else lineageLCA.
			int fastL=(c.level>=0 ? c.level : strL);
			total++;
			if(fastL==treeL){fastMatch++;}
			if(treeL==strL){match++; continue;}
			String key=strL+"->"+treeL;
			mism.merge(key, 1L, Long::sum);
			ArrayList<String> l=samp.get(key);
			if(l==null){samp.put(key, l=new ArrayList<String>());}
			if(l.size()<8){l.add("  tid "+c.taxID+" "+c.name+"\n    "+lin);}
		}
		System.out.println("===== SELF-LCA over ALL records (fresh lineage lineageLCA vs commonAncestorLevel) =====");
		System.out.println("records="+total+"  MATCH(lineage-only)="+match+" = "+fmt(match/(double)total)
			+"   mismatch="+(total-match));
		System.out.println("records="+total+"  MATCH(runtime fast-path)="+fastMatch+" = "+fmt(fastMatch/(double)total)
			+"   mismatch="+(total-fastMatch)+"   <- the ACTUAL sketchLCA self path");
		ArrayList<String> keys=new ArrayList<String>(mism.keySet());
		keys.sort((x,y)->Long.compare(mism.get(y), mism.get(x)));
		for(String key : keys){
			String[] st=key.split("->");
			int sl=Integer.parseInt(st[0]), tl=Integer.parseInt(st[1]);
			System.out.println("\n["+lvlName(sl)+" (str) -> "+lvlName(tl)+" (tree)] count="+mism.get(key));
			for(String s : samp.get(key)){System.out.println(s);}
		}
	}

	/** encodeLevel9 EXACTLY as CladeRanking: invalid (level<0) -> 0; else clamp((11-level)/9). */
	private static float enc(int level){
		if(level<0){return 0f;}
		float x=(11-level)/9f;
		return x<0 ? 0 : (x>1 ? 1 : x);
	}

	private static int idx(int level){return level<0 ? 0 : Math.min(12, level+1);}

	private static boolean colUsed(long[][] m, int col){
		for(int r=0; r<m.length; r++){if(m[r][col]>0){return true;}}
		return false;
	}

	/** level index (matrix index-1): -1=invalid, 0=NO_RANK,1=SUBSP,2=SPECIES...10=DOMAIN,11=LIFE. */
	private static String lvlName(int level){
		switch(level){
			case -1: return "NONE(-1)";
			case 0: return "NORANK(0)";
			case 1: return "SUBSP(1)";
			case 2: return "SPECIES(2)";
			case 3: return "GENUS(3)";
			case 4: return "FAMILY(4)";
			case 5: return "ORDER(5)";
			case 6: return "CLASS(6)";
			case 7: return "PHYLUM(7)";
			case 8: return "KINGDOM(8)";
			case 9: return "SUPERK(9)";
			case 10: return "DOMAIN(10)";
			case 11: return "LIFE(11)";
			default: return "?("+level+")";
		}
	}

	private static String token(String lineage, String prefix){
		String s=";"+lineage;
		int pos=s.indexOf(prefix);
		if(pos<0){return null;}
		int start=pos+prefix.length();
		int end=s.indexOf(';', start);
		return s.substring(start, end<0 ? s.length() : end);
	}

	private static String lname(TaxTree tree, int tid){
		TaxNode tn=tree.getNode(tid, true);
		return tn==null ? "?" : tn.name+" @"+lvlName(tn.level);
	}

	private static String fmt(double d){return String.format("%.5f", d);}
}
