package clade;

import java.io.PrintStream;
import java.util.ArrayList;

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
import shared.Tools;
import structures.ByteBuilder;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Builds neural-network training vectors for QuickClade confidence from a QuickClade
 * machine-format hit TSV. One 48-dimension input vector per hit (main 26 + two 11-dim
 * alt-hit vectors), plus a label. A vector is emitted only for the top MAX_EMIT hits of a
 * query; the alt-hit features look across the top TOP_EXAMINE hits.
 *   mode=continuous: one file; label = LCA(queryTruth, hit) as a 0-1 depth (last column).
 *   mode=binary:     one file PER taxonomic level (out is a prefix); label = 1 if the hit
 *                    is correct at that level (LCA level <= L), else 0.
 *
 * Query truth taxID is parsed from the shred header (tid_NNN or tid|NNN|). A query whose header
 * carries no parseable truth taxID is SKIPPED (fail closed) rather than given a fabricated label.
 * Columns are read BY NAME from the '#'-header (not by hardcoded index), so a flag-driven column
 * shift does not silently miscompute. See the quickclade_nn skill for the full vector spec.
 *
 * Missingness: kid/wkid are printed as -1 by QuickClade when no DDL comparison ran, so a per-hit
 * "sketch-present" flag distinguishes that from a legitimately-measured identity of 0. Likewise a
 * "spectra-present" flag guards the composition/k-mer difs (which default to 1.0 when uncomputed).
 *
 * @author Noire, Brian Bushnell
 */
public class ConfidenceVectorizer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		ConfidenceVectorizer x=new ConfidenceVectorizer(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public ConfidenceVectorizer(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{
			Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			in1=parser.in1;
			out1=parser.out1;
		}
		assert(in1!=null) : "No input (in=machine.tsv).";
		assert(out1!=null) : "No output (out=vectors.tsv or a prefix for mode=binary).";
		assert(treeFile!=null) : "No taxonomy tree (tree=tree.taxtree.gz).";
		assert(buckets>=0) : "buckets must be >=0 (0 = fast/no-sketch mode; else the DB bucket count, e.g. 4096, 32768).";

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
			}else if(a.equals("mode")){
				continuous=(b!=null && (b.startsWith("cont")||b.equals("c")||b.equals("regression")||b.equals("reg")));
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
		outstream.println("featureVersion=2 mode="+(continuous?"continuous":"binary")+" buckets="+buckets+" dims="+VECTOR_DIMS);

		//Set up output writer(s)
		if(continuous){
			writers=new ByteStreamWriter[]{makeBSW(out1)};
		}else{
			writers=new ByteStreamWriter[LEVELS.length];
			for(int i=0; i<LEVELS.length; i++){
				writers[i]=makeBSW(insertSuffix(out1, "_"+LEVEL_NAMES[i]));
			}
		}
		//#dims header: <inputDims> <outputDims>
		//ml.DataLoader parses the #dims header on TAB (DataLoader.delimiter='\t'); space-delimited crashes it.
		ByteBuilder hdr=new ByteBuilder().append("#dims").tab().append(VECTOR_DIMS).tab().append(1).nl();
		for(ByteStreamWriter w : writers){w.print(hdr.toBytes());}

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		LineParser1 lp=new LineParser1('\t');
		byte[] line=bf.nextLine();

		//Header line -> column-name map
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
				flushQuery(group);
				group.clear();
			}
			curQuery=h.qname;
			group.add(h);
		}
		flushQuery(group);

		errorState|=bf.close();
		for(ByteStreamWriter w : writers){if(w!=null){errorState|=w.poisonAndWait();}}
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
		}
		assert(cQName>=0 && cQBases>=0 && cRTid>=0 && cGC>=0 && cK5>=0)
			: "Machine header missing required columns.";
	}

	private Hit parseHit(LineParser1 lp){
		if(lp.terms()<=cRTid){return null;}
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
		return h;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Per-query vectors     ----------------*/
	/*--------------------------------------------------------------*/

	private void flushQuery(ArrayList<Hit> hits){
		if(hits.isEmpty()){return;}
		queries++;
		hitsIn+=hits.size();
		//Fail closed: a shred with no parseable truth taxID (tid_NNN) must NOT get a fabricated label.
		int qTid=BinObject.parseTaxID(hits.get(0).qname);
		if(qTid<1){noLabelQueries++; return;}
		//Emit a vector only for the top MAX_EMIT hits (the strongest candidates, already confidence-sorted).
		int emit=Math.min(MAX_EMIT, hits.size());
		for(int i=0; i<emit; i++){
			float[] v=buildVector(i, hits);
			int caLevel=formalLevel(tree.commonAncestor(qTid, hits.get(i).rtid));
			if(continuous){
				writeRow(writers[0], v, encodeLevel(caLevel));
			}else{
				for(int li=0; li<LEVELS.length; li++){
					float label=(caLevel>0 && caLevel<=LEVELS[li]) ? 1f : 0f;
					writeRow(writers[li], v, label);
				}
			}
			vectorsOut++;
		}
	}

	/** Assemble the 48-dim input vector for hit i. */
	private float[] buildVector(int i, ArrayList<Hit> hits){
		Hit h=hits.get(i);
		float[] v=new float[VECTOR_DIMS];
		int p=0;
		//Per-hit validity of the two composition/sketch methods.
		boolean spectraValid=spectraPresent(h);
		boolean sketchValid=(h.kid>=0); //machine output prints kid/wkid = -1 when no DDL comparison ran
		//--- main (26) ---
		float len=Math.max(h.qBases, 1);
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
		boolean ssuValid=(h.ssuID>0);
		v[p++]=ssuValid ? h.ssuID : 0;                  //15 SSU-ANI
		v[p++]=ssuValid ? 1 : 0;                        //16 SSU-valid
		p=domainOneHot(h.rtid, v, p);                   //17-23 domain 1-hot (7)
		boolean lcaValid=(h.sketchLCALevel>0);
		v[p++]=lcaValid ? encodeLevel(h.sketchLCALevel) : 0; //24 cross-method LCA
		v[p++]=lcaValid ? 1 : 0;                        //25 cross-method LCA valid
		v[p++]=log2buckets;                             //26 log2(#buckets)/16
		//--- small #1: top hit other than this (strongest competitor) ---
		p=smallVector(altTopOther(i, hits), h, hits, v, p);
		//--- small #2: most taxonomically remote hit among the top TOP_EXAMINE ---
		p=smallVector(altMostRemote(i, hits), h, hits, v, p);
		assert(p==VECTOR_DIMS) : p+" != "+VECTOR_DIMS;
		return v;
	}

	/** 11-dim alt-hit vector (or all-invalid if alt<0). Returns new offset. */
	private int smallVector(int alt, Hit self, ArrayList<Hit> hits, float[] v, int p){
		if(alt<0){
			//valid=0; the rest 0
			for(int k=0; k<SMALL_DIMS; k++){v[p+k]=0;}
			return p+SMALL_DIMS;
		}
		Hit a=hits.get(alt);
		boolean aSketch=(a.kid>=0);
		v[p++]=1;                                       //1 valid (alt exists)
		v[p++]=a.k4dif;                                 //2 k4dif
		v[p++]=aSketch ? a.kid : 0;                     //3 kid
		v[p++]=aSketch ? a.wkid : 0;                    //4 wkid
		v[p++]=aSketch ? 1 : 0;                         //5 sketch-present(alt)
		int lcaSelf=formalLevel(tree.commonAncestor(a.rtid, self.rtid));
		v[p++]=encodeLevel(lcaSelf);                    //6 LCA(alt, this) -- always valid
		boolean altLcaValid=(a.sketchLCALevel>0);
		v[p++]=altLcaValid ? encodeLevel(a.sketchLCALevel) : 0; //7 alt's cross-method LCA
		v[p++]=altLcaValid ? 1 : 0;                     //8 valid for d7
		boolean altSsuValid=(a.ssuID>0);
		v[p++]=altSsuValid ? a.ssuID : 0;               //9 SSU-ANI(alt)
		v[p++]=altSsuValid ? 1 : 0;                     //10 SSU-valid(alt)
		v[p++]=a.gcdif;                                 //11 gcdif(query, alt-ref)
		return p;
	}

	/** Index of the top hit other than i: hit1 if i!=0, else hit2(index 1). -1 if none. */
	private int altTopOther(int i, ArrayList<Hit> hits){
		if(hits.size()<2){return -1;}
		return (i==0) ? 1 : 0;
	}

	/** Index of the hit whose LCA with hit i is shallowest (most remote), among the top TOP_EXAMINE
	 *  hits; ties -> highest rank (lowest index, which is strongest since hits are confidence-sorted). */
	private int altMostRemote(int i, ArrayList<Hit> hits){
		Hit self=hits.get(i);
		int best=-1, bestLevel=-1;
		int examine=Math.min(TOP_EXAMINE, hits.size());
		for(int j=0; j<examine; j++){
			if(j==i){continue;}
			int lvl=formalLevel(tree.commonAncestor(self.rtid, hits.get(j).rtid));
			if(lvl>bestLevel){bestLevel=lvl; best=j;} //higher level = shallower/more remote LCA
		}
		return best;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Taxonomy helpers      ----------------*/
	/*--------------------------------------------------------------*/

	/** Spectra comparison ran for this hit? The k-mer difs default to 1.0 (uncomputed); a real
	 *  comparison drives at least one below 1. Rare to be false (spectra is QuickClade's primary method). */
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

	/** Continuous 0-1 LCA encoding: life(11)->0 ... species(2)->1.0, constant step 1/9. */
	private static float encodeLevel(int level){
		float x=(11-level)/9f;
		return x<0 ? 0 : (x>1 ? 1 : x);
	}

	/** Parse a rank-name string ("species".."domain", or ".") to a TaxTree level; -1 if none. */
	private static int levelFromString(String s){
		if(s==null || s.length()==0 || s.equals(".")){return -1;}
		try{return TaxTree.stringToLevel(s.toLowerCase());}
		catch(Throwable e){return -1;}
	}

	/** 7-way domain one-hot of the reference hit: bacteria, archaea, virus, animal, plant, fungi, other-euk. */
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
			if(kn.contains("metazoa")){v[base+3]=1;}          //animal
			else if(kn.contains("viridiplantae")){v[base+4]=1;} //plant
			else if(kn.contains("fungi")){v[base+5]=1;}        //fungi
			else{v[base+6]=1;}                                 //other-eukaryote
		}//else: all zero (unclassified superkingdom)
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

	/** Insert a suffix before the (possibly-compressed) extension: a.tsv.gz + _species -> a_species.tsv.gz */
	private static String insertSuffix(String path, String suffix){
		String ext="";
		String base=path;
		if(base.endsWith(".gz")){ext=".gz"; base=base.substring(0, base.length()-3);}
		int dot=base.lastIndexOf('.');
		if(dot>=0){return base.substring(0,dot)+suffix+base.substring(dot)+ext;}
		return base+suffix+ext;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Hit              ----------------*/
	/*--------------------------------------------------------------*/

	private static class Hit {
		String qname;
		long qBases;
		int rtid;
		float gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif;
		float ssuID=-1, wkid=-1, kid=-1;
		int sketchLCALevel=-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null, out1=null, treeFile=null;
	private int buckets=0;
	private boolean continuous=false;
	private float log2buckets=0;

	private TaxTree tree;
	private ByteStreamWriter[] writers;
	private final FileFormat ffin1;

	private boolean parsedHeader=false;
	//Column indices (by header name); -1 = absent.
	private int cQName=-1, cQBases=-1, cRTid=-1, cGC=-1, cSTR=-1, cHH=-1, cCAGA=-1,
			cK3=-1, cK4=-1, cK5=-1, cSSU=-1, cWKID=-1, cKID=-1, cSketchLCA=-1;

	private long queries=0, hitsIn=0, vectorsOut=0, badRows=0, noLabelQueries=0;

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true, append=false;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final double LN2=Math.log(2);
	static final int MAIN_DIMS=26, SMALL_DIMS=11, VECTOR_DIMS=MAIN_DIMS+2*SMALL_DIMS; //48
	static final int MAX_EMIT=5;    //emit a vector only for the top-5 hits per query
	static final int TOP_EXAMINE=10; //alt-hit features look across the top-10 hits
	//The 9 confidence levels (species..domain). Includes SUPERKINGDOM (=9): the reliable high rank
	//(bacteria/archaea/eukaryota/viruses); NCBI 'kingdom' is sparse/deprecated but kept for completeness.
	static final int[] LEVELS={TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY, TaxTree.ORDER,
			TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.SUPERKINGDOM, TaxTree.DOMAIN};
	static final String[] LEVEL_NAMES={"species","genus","family","order","class","phylum",
			"kingdom","superkingdom","domain"};

}
