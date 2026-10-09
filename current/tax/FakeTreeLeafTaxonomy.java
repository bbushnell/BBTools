package tax;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Assigns real taxonomy parents for permanent FakeTree leaf nodes.
 * @author Brian Bushnell, Yelan
 * @date October 7, 2026
 */
public class FakeTreeLeafTaxonomy {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the FakeTree leaf taxonomy assignment workflow.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		FakeTreeLeafTaxonomy x=new FakeTreeLeafTaxonomy(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public FakeTreeLeafTaxonomy(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("leafmap") || a.equals("map")){
				leafMapFile=b;
			}else if(a.equals("ssu")){
				ssuFile=b;
			}else if(a.equals("quickclade") || a.equals("qc")){
				quickCladeFile=b;
			}else if(a.equals("out")){
				outFile=b;
			}else if(a.equals("outinspect") || a.equals("inspect")){
				inspectFile=b;
			}else if(a.equals("outnodes") || a.equals("nodes")){
				nodeFile=b;
			}else if(a.equals("summary")){
				summaryFile=b;
			}else if(a.equals("tree") || a.equals("treefile") || a.equals("taxtree")){
				treeFile=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				throw new IllegalArgumentException("Unknown parameter "+arg);
			}
		}

		overwrite=parser.overwrite;
		append=parser.append;

		validateParams();

		ffLeafMap=FileFormat.testInput(leafMapFile, FileFormat.TEXT, null, true, true);
		ffSSU=FileFormat.testInput(ssuFile, FileFormat.TEXT, null, true, true);
		ffQuickClade=FileFormat.testInput(quickCladeFile, FileFormat.TEXT, null, true, true);
		ffOut=FileFormat.testOutput(outFile, FileFormat.TEXT, null, true, overwrite, append, false);
		ffInspect=FileFormat.testOutput(inspectFile, FileFormat.TEXT, null, true, overwrite, append, false);
		ffNodes=FileFormat.testOutput(nodeFile, FileFormat.TEXT, null, true, overwrite, append, false);
		ffSummary=FileFormat.testOutput(summaryFile, FileFormat.TEXT, null, true, overwrite, append, false);
	}

	private void validateParams(){
		if(leafMapFile==null){throw new RuntimeException("A leaf map is required: in=<file>");}
		if(ssuFile==null){throw new RuntimeException("An SSU top-hit table is required: ssu=<file>");}
		if(quickCladeFile==null){throw new RuntimeException("A QuickClade top-hit table is required: quickclade=<file>");}
		if(outFile==null){throw new RuntimeException("A proposal output is required: out=<file>");}
		if(nodeFile==null){throw new RuntimeException("A node output is required: outnodes=<file>");}
		if(!Tools.testInputFiles(false, true, leafMapFile, ssuFile, quickCladeFile, treeFile)){
			throw new RuntimeException("Can't read at least one input file.");
		}
		if(!Tools.testOutputFiles(overwrite, append, false, outFile, inspectFile, nodeFile, summaryFile)){
			throw new RuntimeException("Can't write at least one output file.");
		}
		if(!Tools.testForDuplicateFiles(true, leafMapFile, ssuFile, quickCladeFile, treeFile,
				outFile, inspectFile, nodeFile, summaryFile)){
			throw new RuntimeException("Duplicate input/output files are not allowed.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		if(treeFile!=null){TaxTree.treeFile=treeFile;}
		tree=TaxTree.loadTaxTree(treeFile, outstream, false, false);
		if(tree==null){throw new RuntimeException("A TaxTree is required.");}

		ArrayList<Row> leaves=readTable(ffLeafMap);
		HashMap<String, Row> ssu=readTableMap(ffSSU);
		HashMap<String, Row> quickClade=readTableMap(ffQuickClade);

		Collections.sort(leaves, Row.FAKE_TAXID_COMPARATOR);
		ArrayList<Row> proposals=new ArrayList<Row>(leaves.size());
		ArrayList<Row> nodeRows=new ArrayList<Row>(leaves.size());
		ArrayList<Row> inspectRows=new ArrayList<Row>();

		HashSet<Integer> fakeTaxids=new HashSet<Integer>();
		for(Row leaf : leaves){
			String fakeTaxid=leaf.get("fake_taxid");
			boolean added=fakeTaxids.add(Integer.parseInt(fakeTaxid));
			assert(added) : "Duplicate fake_taxid in leaf map: "+fakeTaxid;
			Row row=makeRow(leaf, ssu.get(fakeTaxid), quickClade.get(fakeTaxid));
			proposals.add(row);
			Row nodeRow=makeNodeRow(leaf, row);
			nodeRows.add(nodeRow);
			if(row.get("decision_class").endsWith("conflict") || domainOrHigher(row.get("proposed_parent_rank"))){
				inspectRows.add(row);
			}
			increment(row.get("decision_class"));
			increment("source_"+row.get("evidence_source"));
			increment("node_rank_"+row.get("fake_node_rank"));
		}

		writeRows(proposals, PROPOSAL_FIELDS, ffOut);
		if(ffInspect!=null){writeRows(inspectRows, PROPOSAL_FIELDS, ffInspect);}
		writeRows(nodeRows, NODE_FIELDS, ffNodes);
		if(ffSummary!=null){writeSummary(leaves.size(), inspectRows.size(), ffSummary);}

		t.stop();
		outstream.println("Leaf rows:\t"+leaves.size());
		outstream.println("Rows needing inspection:\t"+inspectRows.size());
		outstream.println("Time:\t"+t);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	private Row makeRow(final Row leaf, final Row srow, final Row qrow){
		final Integer ssuTid=taxid(srow==null ? "" : srow.get("hit_tid"));
		final Integer qcTid=taxid(qrow==null ? "" : qrow.get("ref_taxid"));
		final double ssuAni=(srow==null ? -1 : Double.parseDouble(srow.get("ani")));
		final String ssuTarget=(srow==null ? "" : ssuRank(ssuAni));
		final Parent ssuParent=(ssuTid==null ? Parent.EMPTY : parentAt(ssuTid, ssuTarget));
		final String qcTarget=(qrow==null ? "" : qcRank(qrow));
		final Parent qcParent=(qcTid==null ? Parent.EMPTY : parentAt(qcTid, qcTarget));
		final Parent agreement=deepestAgreement(ssuTid, qcTid);
		final Parent proposed;
		final String source, decision;
		if(srow!=null && ssuAni>=0.800){
			proposed=ssuParent;
			source="findssu";
			decision=(agreement.tid>0 && !"domain".equals(agreement.rank) ? 
				"ssu_with_qc_real_lineage_support" : qrow!=null ? "ssu_with_qc_conflict" : "ssu_only");
		}else if(qrow!=null && qcParent.tid!=TaxTree.LIFE_ID){
			proposed=qcParent;
			source="quickclade";
			decision="quickclade_only";
		}else{
			proposed=fallbackParent(leaf.get("domain"));
			source="fallback";
			decision="unresolved";
		}
		assert(proposed.tid>0) : "No proposed parent for fake_taxid "+leaf.get("fake_taxid")+
			"; SSU="+(srow==null ? "null" : srow.get("hit_tid"))+", QC="+
			(qrow==null ? "null" : qrow.get("ref_taxid"));

		String expectedPhylum=leaf.get("phylum");
		String ssuHitPhylum=srow==null ? "" : srow.get("hit_phylum");
		String qcHitPhylum=qrow==null ? "" : qrow.get("qc_hit_phylum");
		Row out=new Row();
		out.put("fake_taxid", leaf.get("fake_taxid"));
		out.put("collection", leaf.get("collection"));
		out.put("origin", leaf.get("origin"));
		out.put("source_id", leaf.get("source_id"));
		out.put("genome_id", leaf.get("genome_id"));
		out.put("existing_parent_taxid", leaf.get("parent_taxid"));
		out.put("expected_domain", leaf.get("domain"));
		out.put("expected_phylum", expectedPhylum);
		out.put("ssu_hit_tid", srow==null ? "" : srow.get("hit_tid"));
		out.put("ssu_hit_name", srow==null ? "" : srow.get("hit_name"));
		out.put("ssu_hit_phylum", ssuHitPhylum);
		out.put("ssu_ani", srow==null ? "" : srow.get("ani"));
		out.put("ssu_wkid", srow==null ? "" : srow.get("wkid"));
		out.put("ssu_matches", srow==null ? "" : srow.get("matches"));
		out.put("ssu_target_rank", ssuTarget);
		out.put("ssu_parent_taxid", ssuParent.tidString());
		out.put("ssu_parent_rank", ssuParent.rank);
		out.put("ssu_parent_name", ssuParent.name);
		out.put("qc_ref_taxid", qrow==null ? "" : qrow.get("ref_taxid"));
		out.put("qc_ref_name", qrow==null ? "" : qrow.get("ref_name"));
		out.put("qc_hit_phylum", qcHitPhylum);
		out.put("qc_sketch_lca", qrow==null ? "" : qrow.get("sketch_lca"));
		out.put("qc_conf_level", qrow==null ? "" : qrow.get("conf_level"));
		out.put("qc_ani", qrow==null ? "" : qrow.get("ANI"));
		out.put("qc_parent_taxid", qcParent.tidString());
		out.put("qc_parent_rank", qcParent.rank);
		out.put("qc_parent_name", qcParent.name);
		out.put("expected_ssu_phylum_exact", boolString(known(expectedPhylum) && expectedPhylum.equals(ssuHitPhylum)));
		out.put("expected_qc_phylum_exact", boolString(known(expectedPhylum) && expectedPhylum.equals(qcHitPhylum)));
		out.put("ssu_qc_phylum_exact", boolString(known(ssuHitPhylum) && ssuHitPhylum.equals(qcHitPhylum)));
		out.put("deepest_real_agree_taxid", agreement.tidString());
		out.put("deepest_real_agree_rank", agreement.rank);
		out.put("deepest_real_agree_name", agreement.name);
		out.put("proposed_parent_taxid", proposed.tidString());
		out.put("proposed_parent_rank", proposed.rank);
		out.put("proposed_parent_name", proposed.name);
		out.put("fake_node_rank", fakeNodeRank(proposed.rank));
		out.put("evidence_source", source);
		out.put("decision_class", decision);
		return out;
	}

	private Row makeNodeRow(final Row leaf, final Row proposal){
		String rank=proposal.get("fake_node_rank");
		Row out=new Row();
		out.put("fake_taxid", leaf.get("fake_taxid"));
		out.put("parent_taxid", proposal.get("proposed_parent_taxid"));
		out.put("rank", rank);
		out.put("level", rank.equals("subspecies") ? TaxTree.SUBSPECIES : TaxTree.SPECIES);
		out.put("level_extended", rank.equals("subspecies") ? TaxTree.SUBSPECIES_E : TaxTree.SPECIES_E);
		out.put("name", leaf.get("collection")+"__"+leaf.get("source_id"));
		out.put("path_key", leaf.get("collection")+"|"+leaf.get("origin")+"|"+leaf.get("source_id"));
		out.put("collection", leaf.get("collection"));
		out.put("origin", leaf.get("origin"));
		out.put("source_id", leaf.get("source_id"));
		out.put("evidence_source", proposal.get("evidence_source"));
		out.put("decision_class", proposal.get("decision_class"));
		return out;
	}

	private Parent parentAt(int tid, String targetRank){
		TaxNode tn=tree.getNode(tid);
		if(tn==null){return Parent.LIFE;}
		int targetIndex=rankIndex(targetRank);
		while(tn!=null){
			int index=rankIndex(tn.levelStringExtended(false));
			if(index>=targetIndex){return new Parent(tn.id, RANKS[index], tn.name);}
			if(tn.id==tn.pid){break;}
			tn=tree.getNode(tn.pid);
		}
		return Parent.LIFE;
	}

	private Parent deepestAgreement(Integer ssuTid, Integer qcTid){
		if(ssuTid==null || qcTid==null){return Parent.EMPTY;}
		TaxNode ssuNode=tree.getNode(ssuTid), qcNode=tree.getNode(qcTid);
		if(ssuNode==null || qcNode==null){return Parent.EMPTY;}
		HashMap<String, Parent> lineage=makeLineage(ssuNode);
		while(qcNode!=null){
			int index=rankIndex(qcNode.levelStringExtended(false));
			if(index>=0){
				Parent ssu=lineage.get(RANKS[index]);
				if(ssu!=null && ssu.tid==qcNode.id){return ssu;}
			}
			if(qcNode.id==qcNode.pid){break;}
			qcNode=tree.getNode(qcNode.pid);
		}
		return Parent.EMPTY;
	}

	private HashMap<String, Parent> makeLineage(TaxNode tn){
		HashMap<String, Parent> map=new HashMap<String, Parent>();
		while(tn!=null){
			int index=rankIndex(tn.levelStringExtended(false));
			if(index>=0 && !map.containsKey(RANKS[index])){
				map.put(RANKS[index], new Parent(tn.id, RANKS[index], tn.name));
			}
			if(tn.id==tn.pid){break;}
			tn=tree.getNode(tn.pid);
		}
		return map;
	}

	private static ArrayList<Row> readTable(FileFormat ff){
		ArrayList<Row> rows=new ArrayList<Row>();
		ByteFile bf=ByteFile.makeByteFile(ff);
		byte[] line=bf.nextLine();
		if(line==null){throw new RuntimeException("Empty table: "+ff.name());}
		String[] header=new String(line).split("\t", -1);
		for(line=bf.nextLine(); line!=null; line=bf.nextLine()){
			String[] split=new String(line).split("\t", -1);
			Row row=new Row();
			for(int i=0; i<header.length; i++){row.put(header[i], i<split.length ? split[i] : "");}
			rows.add(row);
		}
		boolean error=bf.close();
		if(error){throw new RuntimeException("Error reading "+ff.name());}
		return rows;
	}

	private static HashMap<String, Row> readTableMap(FileFormat ff){
		ArrayList<Row> rows=readTable(ff);
		HashMap<String, Row> map=new HashMap<String, Row>();
		for(Row row : rows){
			map.put(row.get("fake_taxid"), row);
		}
		return map;
	}

	private static void writeRows(ArrayList<Row> rows, String[] fields, FileFormat ff){
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		ByteBuilder bb=new ByteBuilder(1<<16);
		appendHeader(bb, fields);
		for(Row row : rows){
			for(int i=0; i<fields.length; i++){
				if(i>0){bb.tab();}
				bb.append(row.get(fields[i]));
			}
			bb.nl();
			if(bb.length()>=(1<<20)){
				bsw.print(bb);
				bb.clear();
			}
		}
		if(bb.length()>0){bsw.print(bb);}
		boolean error=bsw.poisonAndWait();
		if(error){throw new RuntimeException("Error writing "+ff.name());}
	}

	private void writeSummary(int rows, int inspectRows, FileFormat ff){
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		ByteBuilder bb=new ByteBuilder();
		bb.append("FakeTree leaf-only taxonomy proposal summary\n");
		bb.append("============================================\n\n");
		bb.append("Leaf rows: ").append(rows).nl();
		bb.append("Rows needing inspection: ").append(inspectRows).append("\n\n");
		for(String key : counters.keySet()){
			bb.append(key).tab().append(counters.get(key).longValue()).nl();
		}
		bsw.print(bb);
		boolean error=bsw.poisonAndWait();
		if(error){throw new RuntimeException("Error writing "+ff.name());}
	}

	private static void appendHeader(ByteBuilder bb, String[] fields){
		for(int i=0; i<fields.length; i++){
			if(i>0){bb.tab();}
			bb.append(fields[i]);
		}
		bb.nl();
	}

	private void increment(String key){
		Long old=counters.get(key);
		counters.put(key, old==null ? 1L : old+1);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean domainOrHigher(String rank){return "life".equals(rank) || "domain".equals(rank);}

	private static boolean known(String s){
		return s!=null && s.length()>0 && !s.equals("NA") && !s.equals("None") && !s.equals(".") &&
			!s.equals("-1") && !s.equals("unknown");
	}

	private static String boolString(boolean b){return b ? "true" : "false";}

	private static Integer taxid(String value){
		if(value==null || value.length()<1){return null;}
		for(int i=0; i<value.length(); i++){
			char c=value.charAt(i);
			if(c<'0' || c>'9'){return null;}
		}
		int tid=Integer.parseInt(value);
		return tid>0 ? tid : null;
	}

	private static String ssuRank(double ani){
		if(ani>=0.985){return "genus";}
		if(ani>=0.950){return "family";}
		if(ani>=0.900){return "order";}
		if(ani>=0.800){return "phylum";}
		return "domain";
	}

	private static String qcRank(Row row){
		String lca=row.get("sketch_lca");
		if(rankIndex(lca)>=0){return lca;}
		String conf=row.get("conf_level");
		int colon=conf.indexOf(':');
		if(colon>0){
			String rank=conf.substring(0, colon);
			if(rankIndex(rank)>=0){return rank;}
		}
		return "phylum";
	}

	private static String fakeNodeRank(String parentRank){return "species".equals(parentRank) ? "subspecies" : "species";}

	private static int rankIndex(String rank){
		for(int i=0; i<RANKS.length; i++){
			if(RANKS[i].equals(rank)){return i;}
		}
		return -1;
	}

	private static Parent fallbackParent(String domain){
		if("Bacteria".equals(domain)){return BACTERIA;}
		if("Archaea".equals(domain)){return ARCHAEA;}
		return Parent.LIFE;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static class Row{

		String get(String key){
			String s=map.get(key);
			return s==null ? "" : s;
		}

		void put(String key, Object value){map.put(key, value==null ? "" : value.toString());}

		private final HashMap<String, String> map=new HashMap<String, String>();

		static final Comparator<Row> FAKE_TAXID_COMPARATOR=new Comparator<Row>(){
			@Override
			public int compare(Row a, Row b){
				return Integer.parseInt(a.get("fake_taxid"))-Integer.parseInt(b.get("fake_taxid"));
			}
		};

	}

	private static class Parent{

		Parent(int tid_, String rank_, String name_){
			tid=tid_;
			rank=rank_;
			name=name_;
		}

		String tidString(){return tid<1 ? "" : Integer.toString(tid);}

		final int tid;
		final String rank, name;

		static final Parent EMPTY=new Parent(-1, "", "");
		static final Parent LIFE=new Parent(TaxTree.LIFE_ID, "life", "Life");

	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String leafMapFile=null;
	private String ssuFile=null;
	private String quickCladeFile=null;
	private String outFile=null;
	private String inspectFile=null;
	private String nodeFile=null;
	private String summaryFile=null;
	private String treeFile=TaxTree.defaultTreeFile();

	private TaxTree tree;

	private FileFormat ffLeafMap, ffSSU, ffQuickClade;
	private FileFormat ffOut, ffInspect, ffNodes, ffSummary;

	private boolean overwrite=true;
	private boolean append=false;
	private boolean errorState=false;

	private final TreeMap<String, Long> counters=new TreeMap<String, Long>();

	private PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] RANKS={
		"species", "genus", "family", "order", "class", "phylum", "domain"
	};

	private static final String[] PROPOSAL_FIELDS={
		"fake_taxid", "collection", "origin", "source_id", "genome_id",
		"existing_parent_taxid", "expected_domain", "expected_phylum",
		"ssu_hit_tid", "ssu_hit_name", "ssu_hit_phylum", "ssu_ani", "ssu_wkid",
		"ssu_matches", "ssu_target_rank", "ssu_parent_taxid", "ssu_parent_rank",
		"ssu_parent_name", "qc_ref_taxid", "qc_ref_name", "qc_hit_phylum",
		"qc_sketch_lca", "qc_conf_level", "qc_ani", "qc_parent_taxid",
		"qc_parent_rank", "qc_parent_name", "expected_ssu_phylum_exact",
		"expected_qc_phylum_exact", "ssu_qc_phylum_exact",
		"deepest_real_agree_taxid", "deepest_real_agree_rank",
		"deepest_real_agree_name", "proposed_parent_taxid",
		"proposed_parent_rank", "proposed_parent_name", "fake_node_rank",
		"evidence_source", "decision_class"
	};

	private static final String[] NODE_FIELDS={
		"fake_taxid", "parent_taxid", "rank", "level", "level_extended",
		"name", "path_key", "collection", "origin", "source_id",
		"evidence_source", "decision_class"
	};

	private static final Parent BACTERIA=new Parent(TaxTree.BACTERIA_ID, "domain", "Bacteria");
	private static final Parent ARCHAEA=new Parent(TaxTree.ARCHAEA_ID, "domain", "Archaea");

	public static boolean verbose=false;

}
