package clade;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import bin.BinObject;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Builds neural-network training vectors for QuickClade RANKING from a QuickClade machine-format
 * hit TSV. One 33-dimension input vector per hit plus a regression label; the trained net reorders
 * a query's hits, replacing the compositeScore heuristic as the display ordering.
 *
 * Input MUST come from a run with printcomposite=t (the trailing CompositeScore column) and
 * records=10 heapsize=10; the tool fails loudly rather than silently scoring a file without it.
 *
 * VECTOR (33):
 *   1-26   the MAIN block of the 48-dim confidence vector, unchanged. Deliberately DUPLICATED from
 *          ConfidenceVectorizer.buildVector rather than shared, so this tool cannot destabilize a
 *          shipped one; RankingEquivalence is the gate that proves the two have not drifted.
 *   27     (x-t)/((m-t)+EPS) clamped [-2,2] -- composite normalized within the query
 *   28     0.1*(m-t) capped 1.0 -- the spread, i.e. how much dim 27's position should be believed
 *   29     0.1*min(hits in query, 10)
 *   30-33  0.1 * how many of the top TOP_COUNT share species / genus / family / phylum with this hit
 *
 * EPS=0.02 is MEASURED, not chosen for convenience: below a within-query spread of 0.02,
 * compositeScore's ordering is at or below a coin flip in every clade tested (bacteria 0.463,
 * fungi 0.546-0.584), so without EPS the normalization would stretch a pure noise band across the
 * full output range and manufacture an ordering. See plans/ranking_nn.plan section 4.
 *
 * The sharing counts INCLUDE SELF (Brian): self-exclusion makes them sum to 9 sometimes and 10
 * others for no benefit. They replace the confidence vector's two alt-hit blocks, which described
 * a single competitor and swung 11 of 48 dims on one hit; a count of how many of the top hits
 * corroborate this one at each rank is the same information as a density, and degrades gracefully.
 *
 * LABEL: 11-rung LCA ladder, same_node=1.0, species=0.9 ... domain=0.1, life=0. The rung above
 * species exists because strain-level hits all tied at the ceiling under the 10-rung encoding --
 * precisely where ranking decisions are made.
 *
 * Query truth taxID is parsed from the shred header (tid_NNN); a query with no parseable truth is
 * SKIPPED rather than given a fabricated label. Columns are read BY NAME from the '#'-header.
 *
 * @author Noire, Brian Bushnell
 */
public class RankingVectorizer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		RankingVectorizer x=new RankingVectorizer(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public RankingVectorizer(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{
			Parser parser=parse(args);
			in1=parser.in1;
			out1=parser.out1;
		}
		assert(in1!=null) : "No input (in=machine.tsv).";
		assert(out1!=null) : "No output (out=vectors.tsv).";
		assert(treeFile!=null) : "No taxonomy tree (tree=tree.taxtree.gz).";
		assert(buckets>=0) : "buckets must be >=0 (0 = fast/no-sketch mode; else the DB bucket count).";
		assert(eps>=0) : "eps must be >=0.";

		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}

	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
			}else if(a.equals("buckets") || a.equals("ddlbuckets")){
				buckets=Integer.parseInt(b);
			}else if(a.equals("eps")){
				eps=Float.parseFloat(b);
			}else if(a.equals("maxemit")){
				maxEmit=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//handled by parser (in=, out=, overwrite=, ...)
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		outstream.println("Loading tree: "+treeFile);
		tree=TaxTree.loadTaxTree(treeFile, outstream, true, false);
		log2buckets=(float)(Math.log(Math.max(buckets,1))/Math.log(2)/16.0);
		outstream.println("rankingVector dims="+VECTOR_DIMS+" buckets="+buckets+" eps="+eps+
				" maxEmit="+maxEmit+" topCount="+TOP_COUNT);

		ByteStreamWriter bsw=makeBSW(out1);
		//ml.DataLoader parses the #dims header on TAB; space-delimited crashes it.
		bsw.print(new ByteBuilder().append("#dims").tab().append(VECTOR_DIMS).tab().append(1).nl().toBytes());

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		LineParser1 lp=new LineParser1('\t');
		byte[] line=bf.nextLine();

		while(line!=null && (line.length==0 || (line[0]=='#' && !parsedHeader))){
			if(line.length>0 && line[0]=='#'){parseHeader(line, lp); parsedHeader=true; break;}
			line=bf.nextLine();
		}
		assert(parsedHeader) : "No '#'-header line found in "+in1;

		ArrayList<Hit> group=new ArrayList<Hit>();
		String curQuery=null;
		for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			Hit h=parseHit(lp);
			if(h==null){badRows++; continue;}
			if(curQuery!=null && !h.qname.equals(curQuery)){
				flushQuery(group, bsw);
				group.clear();
			}
			curQuery=h.qname;
			group.add(h);
		}
		flushQuery(group, bsw);

		errorState|=bf.close();
		errorState|=bsw.poisonAndWait();
		t.stop();
		outstream.println("Queries: "+queries+"  hits: "+hitsIn+"  vectors: "+vectorsOut+
				"  noLabelQueries: "+noLabelQueries+"  badRows: "+badRows);
		outstream.println("Time: \t"+t);
		if(errorState){throw new RuntimeException(getClass().getName()+" terminated in an error state.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Header / row parsing    ----------------*/
	/*--------------------------------------------------------------*/

	private void parseHeader(byte[] line, LineParser1 lp){
		lp.set(line);
		int n=lp.terms();
		for(int i=0; i<n; i++){
			String name=lp.parseString(i);
			if(name.startsWith("#")){name=name.substring(1);}
			if(name.equals("QueryName")){cQName=i;}
			else if(name.equals("Q_Bases")){cQBases=i;}
			else if(name.equals("R_TaxID")){cRTid=i;}
			else if(name.equals("GCdif")){cGC=i;}
			else if(name.equals("STRdif")){cSTR=i;}
			else if(name.equals("HHdif")){cHH=i;}
			else if(name.equals("CAGAdif")){cCAGA=i;}
			else if(name.equals("k3dif")){cK3=i;}
			else if(name.equals("k4dif")){cK4=i;}
			else if(name.equals("k5dif")){cK5=i;}
			else if(name.equals("ssuID")){cSSU=i;}
			else if(name.equals("WKID")){cWKID=i;}
			else if(name.equals("KID")){cKID=i;}
			else if(name.equals("Sketch_LCA")){cSketchLCA=i;}
			else if(name.equals("CompositeScore")){cComposite=i;}
		}
		assert(cQName>=0 && cQBases>=0 && cRTid>=0 && cGC>=0 && cK5>=0)
			: "Machine header missing required columns.";
		//Crash loud rather than train a ranking net on a file that never carried the ranking signal.
		if(cComposite<0){
			throw new RuntimeException("No CompositeScore column in "+in1+
					" -- regenerate hits with quickclade printcomposite=t");
		}
	}

	private Hit parseHit(LineParser1 lp){
		if(lp.terms()<=cRTid || lp.terms()<=cComposite){return null;}
		Hit h=new Hit();
		h.qname=lp.parseString(cQName);
		h.qBases=lp.parseLong(cQBases);
		h.rtid=lp.parseInt(cRTid);
		if(h.rtid<1){return null;}
		h.gcdif=lp.parseFloat(cGC);
		h.strdif=lp.parseFloat(cSTR);
		h.hhdif=lp.parseFloat(cHH);
		h.cagadif=lp.parseFloat(cCAGA);
		h.k3dif=lp.parseFloat(cK3);
		h.k4dif=lp.parseFloat(cK4);
		h.k5dif=lp.parseFloat(cK5);
		h.ssuID=(cSSU>=0 ? lp.parseFloat(cSSU) : -1);
		h.wkid=(cWKID>=0 ? lp.parseFloat(cWKID) : -1);
		h.kid=(cKID>=0 ? lp.parseFloat(cKID) : -1);
		h.sketchLCALevel=(cSketchLCA>=0 ? levelFromString(lp.parseString(cSketchLCA)) : -1);
		h.composite=lp.parseFloat(cComposite);
		return h;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Per-query vectors     ----------------*/
	/*--------------------------------------------------------------*/

	private void flushQuery(ArrayList<Hit> hitsIn0, ByteStreamWriter bsw){
		if(hitsIn0.isEmpty()){return;}
		queries++;
		hitsIn+=hitsIn0.size();
		int qTid=BinObject.parseTaxID(hitsIn0.get(0).qname);
		if(qTid<1){noLabelQueries++; return;}

		//ONE PASS: m, t and the sharing counts are all DEFINED over the compositeScore ordering.
		//The net then reranks from those fixed inputs; it never feeds back. Two sorts total in the
		//whole system -- composite here, NN score at display time -- and never an iteration.
		ArrayList<Hit> hits=new ArrayList<Hit>(hitsIn0);
		Collections.sort(hits, (a, b) -> Float.compare(b.composite, a.composite));

		final int n=hits.size();
		final float m=hits.get(0).composite;
		//10th place, not min: min drifts with the record count, so the same query scored with
		//records=10 vs records=50 would otherwise normalize to different numbers.
		final float t=hits.get(Math.min(TOP_COUNT-1, n-1)).composite;
		final float spread=m-t;

		final int top=Math.min(TOP_COUNT, n);
		final int emit=Math.min(maxEmit, n);
		//Ancestors are needed for every EMITTED hit, not just the top-10 window: with maxemit>10 a
		//hit outside the window still has a well-defined "how many of the top 10 corroborate me at
		//this rank", and it is usually nonzero. Counting only over the window is not self-exclusion
		//(Brian: no self-exclusion) -- a hit inside the window counts itself because it is IN the
		//set being counted; a hit outside simply is not.
		final int[][] anc=new int[Math.max(top, emit)][SHARE_LEVELS.length];
		for(int i=0; i<anc.length; i++){
			for(int L=0; L<SHARE_LEVELS.length; L++){
				anc[i][L]=ancestorId(hits.get(i).rtid, SHARE_LEVELS[L]);
			}
		}

		for(int i=0; i<emit; i++){
			float[] v=buildVector(i, hits, m, t, spread, anc, top);
			writeRow(bsw, v, labelDepth(qTid, hits.get(i).rtid));
			vectorsOut++;
		}
	}

	/** Assemble the 33-dim input vector for hit i. */
	private float[] buildVector(int i, ArrayList<Hit> hits, float m, float t, float spread,
			int[][] anc, int top){
		final Hit h=hits.get(i);
		final float[] v=new float[VECTOR_DIMS];
		int p=0;

		//--- main (26): byte-for-byte the same construction as ConfidenceVectorizer.buildVector
		//dims 1-26. RankingEquivalence asserts they stay identical. ---
		final boolean spectraValid=spectraPresent(h);
		final boolean sketchValid=(h.kid>=0); //machine output prints kid/wkid = -1 when no DDL ran
		final float len=Math.max(h.qBases, 1);
		v[p++]=(float)(Math.log(len)/LN2);              //1 log2(len)
		v[p++]=(float)(0.001*Math.sqrt(len));           //2 0.001*sqrt(len)
		v[p++]=spectraValid ? h.gcdif : 0;              //3 gcdif
		v[p++]=spectraValid ? h.strdif : 0;             //4 strdif
		v[p++]=spectraValid ? h.hhdif : 0;              //5 hhdif
		v[p++]=spectraValid ? h.cagadif : 0;            //6 cagadif
		v[p++]=spectraValid ? h.k3dif : 0;              //7 k3dif
		v[p++]=spectraValid ? h.k4dif : 0;              //8 k4dif
		v[p++]=spectraValid ? h.k5dif : 0;              //9 k5dif
		v[p++]=spectraValid ? h.k3dif/(h.k5dif+0.01f) : 0; //10 k3dif/(k5dif+0.01)
		v[p++]=sketchValid ? h.kid : 0;                 //11 kid
		v[p++]=sketchValid ? h.wkid : 0;                //12 wkid
		v[p++]=spectraValid ? 1 : 0;                    //13 spectra-present
		v[p++]=sketchValid ? 1 : 0;                     //14 sketch-present
		final boolean ssuValid=(h.ssuID>0);
		v[p++]=ssuValid ? h.ssuID : 0;                  //15 SSU-ANI
		v[p++]=ssuValid ? 1 : 0;                        //16 SSU-valid
		p=domainOneHot(h.rtid, v, p);                   //17-23 domain 1-hot (7)
		final boolean lcaValid=(h.sketchLCALevel>0);
		v[p++]=lcaValid ? encodeLevel9(h.sketchLCALevel) : 0; //24 cross-method LCA
		v[p++]=lcaValid ? 1 : 0;                        //25 cross-method LCA valid
		v[p++]=log2buckets;                             //26 log2(#buckets)/16
		assert(p==MAIN_DIMS) : p+" != "+MAIN_DIMS;

		//--- ranking block (7) ---
		final float denom=spread+eps;
		float norm=(denom<=0 ? 0 : (h.composite-t)/denom);
		norm=Math.max(-CLAMP, Math.min(CLAMP, norm));
		v[p++]=norm;                                    //27 normalized composite
		v[p++]=Math.min(1f, 0.1f*spread);               //28 spread magnitude
		v[p++]=0.1f*Math.min(hits.size(), TOP_COUNT);   //29 hit count
		//Counted over the top-TOP_COUNT window, so a hit inside it counts itself (>=1) and a hit
		//emitted beyond it does not appear in its own count and may legitimately score 0.
		for(int L=0; L<SHARE_LEVELS.length; L++){      //30-33 species/genus/family/phylum
			v[p++]=0.1f*sharingCount(i, L, anc, top);
		}
		assert(p==VECTOR_DIMS) : p+" != "+VECTOR_DIMS;
		return v;
	}

	/** How many of the top TOP_COUNT hits share hit i's ancestor at share-level L. */
	private static int sharingCount(int i, int L, int[][] anc, int top){
		final int self=anc[i][L];
		if(self<1){return 0;} //no resolvable ancestor at this rank -- corroboration is undefined, not zero-evidence
		int count=0;
		for(int j=0; j<top; j++){
			if(anc[j][L]==self){count++;}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Taxonomy helpers      ----------------*/
	/*--------------------------------------------------------------*/

	/** TaxID of this reference's ancestor at the given rank, or -1 if none. */
	private int ancestorId(int rtid, int level){
		if(rtid<1){return -1;}
		TaxNode n=tree.getNodeAtLevel(rtid, level);
		return n==null ? -1 : n.id;
	}

	/** The ranking label: 11-rung ladder, with same_node above species. */
	private float labelDepth(int qTid, int rTid){
		if(qTid==rTid){return 1f;} //exact organism -- the rung species could never express
		final int level=formalLevel(tree.commonAncestor(qTid, rTid));
		final float x=(11-level)/10f;
		return x<0 ? 0 : (x>1 ? 1 : x);
	}

	private static boolean spectraPresent(Hit h){
		return h.k3dif<1 || h.k4dif<1 || h.k5dif<1;
	}

	/** Promote a node's level up the parent chain to the next formal rank (>=SPECIES). */
	private int formalLevel(int tid){
		if(tid<1){return TaxTree.LIFE;}
		TaxNode n=tree.getNode(tid);
		while(n!=null && n.level<TaxTree.SPECIES){n=(n.pid>0 && n.pid!=n.id) ? tree.getNode(n.pid) : null;}
		return n==null ? TaxTree.LIFE : n.level;
	}

	/** The 10-rung encoding used INSIDE the main block (dim 24), unchanged from the confidence
	 *  vector. Distinct from labelDepth's 11 rungs, which applies only to the label. */
	private static float encodeLevel9(int level){
		float x=(11-level)/9f;
		return x<0 ? 0 : (x>1 ? 1 : x);
	}

	private static int levelFromString(String s){
		if(s==null || s.length()==0 || s.equals(".")){return -1;}
		try{return TaxTree.stringToLevel(s.toLowerCase());}
		catch(Throwable e){return -1;}
	}

	/** 7-way domain one-hot of the reference hit. */
	private int domainOneHot(int rtid, float[] v, int p){
		int base=p;
		for(int k=0; k<7; k++){v[base+k]=0;}
		TaxNode sk=tree.getNodeAtLevel(rtid, TaxTree.SUPERKINGDOM);
		String skn=(sk!=null && sk.name!=null) ? sk.name.toLowerCase() : "";
		if(skn.contains("bacteria")){v[base+0]=1;}
		else if(skn.contains("archaea")){v[base+1]=1;}
		else if(skn.contains("virus") || skn.contains("viroid") || skn.contains("viria")){v[base+2]=1;}
		else if(skn.contains("eukaryota")){
			TaxNode k=tree.getNodeAtLevel(rtid, TaxTree.KINGDOM);
			String kn=(k!=null && k.name!=null) ? k.name.toLowerCase() : "";
			if(kn.contains("metazoa")){v[base+3]=1;}
			else if(kn.contains("viridiplantae")){v[base+4]=1;}
			else if(kn.contains("fungi")){v[base+5]=1;}
			else{v[base+6]=1;}
		}
		return p+7;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeRow(ByteStreamWriter w, float[] v, float label){
		ByteBuilder bb=new ByteBuilder(v.length*10);
		for(int i=0; i<v.length; i++){bb.append(v[i], 7).tab();}
		bb.append(label, 7).nl();
		w.print(bb.toBytes());
	}

	private static ByteStreamWriter makeBSW(String path){
		FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, true, true, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Hit              ----------------*/
	/*--------------------------------------------------------------*/

	static class Hit {
		String qname;
		long qBases;
		int rtid;
		float gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif;
		float ssuID=-1, wkid=-1, kid=-1;
		int sketchLCALevel=-1;
		float composite;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null, out1=null, treeFile=null;
	private int buckets=0;
	private int maxEmit=MAX_EMIT;
	private float eps=DEFAULT_EPS;
	private float log2buckets=0;

	private TaxTree tree;
	private final FileFormat ffin1;

	private boolean parsedHeader=false;
	private int cQName=-1, cQBases=-1, cRTid=-1, cGC=-1, cSTR=-1, cHH=-1, cCAGA=-1,
			cK3=-1, cK4=-1, cK5=-1, cSSU=-1, cWKID=-1, cKID=-1, cSketchLCA=-1, cComposite=-1;

	private long queries=0, hitsIn=0, vectorsOut=0, badRows=0, noLabelQueries=0;

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final double LN2=Math.log(2);
	static final int MAIN_DIMS=26, RANK_DIMS=7, VECTOR_DIMS=MAIN_DIMS+RANK_DIMS; //33
	static final int MAX_EMIT=10;   //emit a vector for the top-10 hits: those are the rows the net must sort
	static final int TOP_COUNT=10;  //m, t and the sharing counts are all defined over the top 10
	static final float CLAMP=2f;    //dim 27 bound; measured to be hit 0.00% of the time, so it is free insurance
	static final float DEFAULT_EPS=0.02f;
	static final int[] SHARE_LEVELS={TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY, TaxTree.PHYLUM};

}
