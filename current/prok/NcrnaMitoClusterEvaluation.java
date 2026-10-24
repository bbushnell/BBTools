package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import tax.TaxNode;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

/** Evaluates a complete native cluster arm without exclusions or merge decisions.
 * Quantum uses the shorter input as query; equal lengths use the first named input.
 * Pair output contains each unordered cluster pair, including diagonals, exactly
 * once in FASTA order. Identities are descriptive, not biological classification.
 * @author Ganyu */
public final class NcrnaMitoClusterEvaluation implements Accumulator<NcrnaMitoClusterEvaluation.Worker> {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoClusterEvaluation.class,false);
		try{final NcrnaMitoClusterEvaluation x=new NcrnaMitoClusterEvaluation();x.parse(pp.args);x.process();
			pp.outstream.println("MITO_CLUSTER_EVALUATION_PASS members="+x.members.size()+" unique="+x.sequences.size()+" models="+x.models.size()+" assigned="+x.assigned+" unassigned="+(x.members.size()-x.assigned)+" alignments="+x.tasks.size()+" no_exclusion=true no_merge_cutoff=true");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();String verified=null;
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("assignments")){assignments=v;}else if(k.equals("consensus")){consensus=v;}else if(k.equals("hbm")){hbm=v;}
			else if(k.equals("audit")){audit=v;}else if(k.equals("tree")){treeFile=v;}else if(k.equals("lineagerank")){rankName=v;}
			else if(k.equals("arm")){arm=v;}else if(k.equals("outprefix")){prefix=v;}else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}
			else if(k.equals("verified")){verified=v;}else if(k.equals("lengthbins")){bins=NcrnaMitoCharacterization.parseBins(v);}
			else if(k.equals("ref16")){refFiles[0]=v;}else if(k.equals("ref18")){refFiles[1]=v;}else if(k.equals("ref23")){refFiles[2]=v;}else if(k.equals("reflsu")){refFiles[3]=v;}
			else if(k.equals("ref16id")){refIds[0]=v;}else if(k.equals("ref18id")){refIds[1]=v;}else if(k.equals("ref23id")){refIds[2]=v;}else if(k.equals("reflsuid")){refIds[3]=v;}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		input=p.in1;require(expected>0 && prefix!=null && arm!=null && rankName!=null && !p.append && p.maxReads<0,"Require an explicit complete arm, expected count, output prefix and lineage rank; no append/read limit");
		require(verified!=null,"Require producer-verified input manifest");final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);
		for(String path:new String[]{input,assignments,consensus,hbm,audit,treeFile}){require(path!=null && pins.get(path)!=null,"Missing pinned input: "+path);}
		for(int i=0;i<4;i++){require(refFiles[i]!=null && refIds[i]!=null && pins.get(refFiles[i])!=null,"All four universal FASTAs and exact record IDs must be explicit and pinned");}
		require(treeFile.indexOf(',')<0,"Supply one explicitly pinned taxonomy tree, not an implicit auxiliary list");
		rank=TaxTree.stringToLevelExtended(rankName);require(rank>0,"Lineage rank must be a named taxonomy rank");
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
		RrnaResourceIO.fresh(prefix,SUFFIXES);
	}
	private void process(){
		RrnaResourceIO.stream(input,r->{final Member m=new Member();m.header=r.id;m.id=RrnaResourceIO.first(r.id);
			require(byId.get(m.id)==null && r.bases!=null && r.length()>0,"Duplicate or empty arm member: "+m.id);
			m.uniqueId=RrnaResourceIO.field(r.id,"unique");m.taxid=TaxTree.parseTaxID(r.id);require(m.taxid>0,"Source member requires positive tid: "+m.id);
			m.partial=flag(RrnaResourceIO.field(r.id,"partial"));m.pseudo=flag(RrnaResourceIO.field(r.id,"pseudo"));
			Sequence s=byUnique.get(m.uniqueId);if(s==null){s=new Sequence(m.uniqueId,r.bases.clone(),sequences.size());byUnique.put(m.uniqueId,s);sequences.add(s);}
			else{require(Arrays.equals(s.bases,r.bases),"Unique source ID binds conflicting bases: "+m.uniqueId);}
			m.sequence=s;members.add(m);byId.put(m.id,m);
		});require(members.size()==expected,"Arm member count differs from its declared complete denominator");
		RrnaResourceIO.stream(consensus,r->{final String token=RrnaResourceIO.first(r.id);require(r.bases!=null && r.length()>0 && modelIds.get(r.id)<0 && modelTokens.get(token)<0,"Empty/duplicate consensus name or first-token model ID: "+r.id);
			modelTokens.put(token,models.size());modelIds.put(r.id,models.size());models.add(new Read(r.bases.clone(),null,r.id,models.size()));});require(!models.isEmpty(),"An evaluated library must contain at least one model");
		readHbm();readAssignments();readAudit();loadLineages();
		for(int i=0;i<4;i++){final int ri=i;RrnaResourceIO.stream(refFiles[i],r->{if(RrnaResourceIO.first(r.id).equals(refIds[ri])){require(refs[ri]==null && r.bases!=null && r.length()>0,"Duplicate/empty selected universal reference");refs[ri]=r.bases.clone();}});require(refs[i]!=null,"Missing exact universal record: "+refIds[i]);}
		for(Member m:members){if(m.model>=0){final String label=models.get(m.model).id;Task t=m.sequence.assigned.get(label);
			if(t==null){t=new Task("MEMBER",m.uniqueId,label,m.sequence.bases,models.get(m.model).bases);m.sequence.assigned.put(label,t);tasks.add(t);}m.task=t;assigned++;}}
		for(int i=0;i<models.size();i++){final Read a=models.get(i);
			for(int u=0;u<4;u++){tasks.add(new Task("UNIVERSAL",a.id,REF_NAMES[u]+":"+refIds[u],a.bases,refs[u]));}
			for(int j=i;j<models.size();j++){final Read b=models.get(j);tasks.add(new Task("CLUSTER_PAIR",a.id,b.id,a.bases,b.bases));}
		}
		final ArrayList<Worker> workers=new ArrayList<Worker>();for(int i=0;i<Math.min(Shared.threads(),tasks.size());i++){workers.add(new Worker(this));}
		require(ThreadWaiter.startAndWait(workers,this),"An alignment worker failed; incomplete measurements cannot become a report");
		write();
	}
	private void readAssignments(){
		final ByteFile in=RrnaResourceIO.open(assignments);final LineParser1 p=new LineParser1('\t');int count=0;final int[] kept=new int[models.size()];
		try{RrnaResourceIO.header(in,"#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus",assignments);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==6,"Native assignment schema changed");final String header=p.parseString(0);final Member m=byId.get(RrnaResourceIO.first(header));
				require(m!=null && !m.assignmentSeen && m.header.equals(header) && m.taxid==p.parseInt(1),"Assignments must conserve exact source headers once");m.assignmentSeen=true;count++;
				m.clusterKey=p.parseString(3);m.status=p.parseString(5);final String label=p.parseString(4);
				if(m.status.equals("kept")){m.model=modelIds.get(label);require(m.model>=0,"Kept assignment has no matching library model: "+label);kept[m.model]++;}
				else{require((m.status.equals("small") || m.status.equals("orphan")) && label.isEmpty(),"Unknown unassigned status or model label: "+m.status);}
			}
		}finally{require(!in.close(),"Assignment read failed");}require(count==expected,"Missing assignment rows must fail, not silently become unassigned");
		for(int i=0;i<models.size();i++){final int declared=Integer.parseInt(RrnaResourceIO.field(models.get(i).id,"n"));require(declared>0 && declared==kept[i],"Declared model n= must equal its actual kept-assignment count: "+models.get(i).id+" versus "+kept[i]);}
	}
	private void readAudit(){
		final ByteFile in=RrnaResourceIO.open(audit);final LineParser1 p=new LineParser1('\t');int found=0;
		try{RrnaResourceIO.header(in,NcrnaMitoCharacterization.HEADER,audit);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==30,"Construction audit schema changed");if(!p.termEquals("ANNOTATION",0)){require(p.termEquals("CM_ONLY_EXACT_GEOMETRY",0),"Unknown audit kind");continue;}
				final Member m=byId.get(p.parseString(1));if(m==null){continue;}
				require(!m.auditSeen && m.sequence.bases.length==p.parseInt(6) && m.taxid==p.parseInt(3) && m.uniqueId.equals(p.parseString(12)) && m.partial==flag(p.parseString(7)) && m.pseudo==flag(p.parseString(8)),"Selected arm and accepted audit disagree for "+m.id);
				m.auditSeen=true;found++;m.product=p.parseString(9);m.preference=p.parseString(15);m.cm=p.parseString(18);m.concordance=p.parseString(23);
				m.bin=NcrnaMitoCharacterization.bin(m.sequence.bases.length,bins);
			}
		}finally{require(!in.close(),"Audit read failed");}require(found==expected,"Every selected member requires exactly one accepted audit row");
	}
	private void readHbm(){
		final ByteFile in=RrnaResourceIO.open(hbm);final LineParser1 p=new LineParser1('\t');int index=0;
		try{RrnaResourceIO.header(in,"#HBM\t1",hbm);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(index<models.size() && p.terms()==2 && p.termStartsWith(">",0) && p.termStartsWith("depth=",1),"Expected one ordered HBM header per consensus");
				final Read model=models.get(index++);require(p.parseString(0).substring(1).equals(model.id),"HBM/FASTA names or order differ");
				final byte[] bases=in.nextLine();require(bases!=null && Arrays.equals(bases,model.bases),"HBM/FASTA bases or coordinate length differ");
				for(int pos=0;pos<bases.length;pos++){final byte[] row=in.nextLine();require(row!=null,"Truncated HBM counts");p.set(row);require(p.terms()==5,"HBM needs A/C/G/T/deletion counts per position");for(int a=0;a<5;a++){require(p.parseInt(a)>=0,"Negative HBM allele count");}}
			}
		}finally{require(!in.close(),"HBM read failed");}require(index==models.size(),"HBM and FASTA model counts differ");
	}
	private void loadLineages(){
		assert(treeFile!=null && rank>0):"Explicit pinned taxonomy and rank precede grouping";
		final boolean old=TaxTree.LOAD_FAKE_TREE;final TaxTree tree;
		try{TaxTree.LOAD_FAKE_TREE=false;tree=TaxTree.loadTaxTreeFile(treeFile);}finally{TaxTree.LOAD_FAKE_TREE=old;}
		final IntHashMap cache=new IntHashMap();final ArrayList<String> labels=new ArrayList<String>();
		for(Member m:members){int i=cache.get(m.taxid);if(i<0){i=labels.size();labels.add(lineage(tree,m.taxid,rank));cache.put(m.taxid,i);}m.lineage=labels.get(i);}
	}
	public static String lineage(TaxTree tree,int tid,int rank){
		require(tree!=null && tid>0 && rank>0,"Taxonomy lookup needs a source taxid and named rank");TaxNode n=tree.getNode(tid,true);if(n==null){return "MISSING_TAXID";}
		final IntList visited=new IntList();
		while(true){require(!visited.contains(n.id),"Cycle in taxonomy parent chain");visited.add(n.id);
			if(n.levelExtended==rank){return n.id+":"+n.name.replace('\t',' ').replace('\n',' ');}
			if(n.id==n.pid){return "RANK_ABSENT";}n=tree.getNode(n.pid,true);if(n==null){return "MISSING_PARENT";}
		}
	}
	private void write(){
		assert(assigned<=expected):"Assigned and explicit unassigned rows share the complete arm denominator";
		final ByteStreamWriter raw=RrnaResourceIO.writer(prefix+SUFFIXES[0]),pairs=RrnaResourceIO.writer(prefix+SUFFIXES[1]),summary=RrnaResourceIO.writer(prefix+SUFFIXES[2]),composition=RrnaResourceIO.writer(prefix+SUFFIXES[3]);
		final ByteBuilder b=new ByteBuilder();final Summary[] overall=summaries();final Summary[][] perModel=new Summary[models.size()+1][];
		final ObjectMap<String,Summary[]> perLineage=new ObjectMap<String,Summary[]>(String.class,Summary[].class);final ArrayList<String> lineageOrder=new ArrayList<String>();
		final Composition[] comp=new Composition[perModel.length];for(int i=0;i<perModel.length;i++){perModel[i]=summaries();comp[i]=new Composition(bins.length+1);}
		try{
			raw.println("arm\tlocus_id\tunique_id\tsource_taxid\tsource_length\tpartial\tpseudo\tsource_product\tlength_bin\tlineage_rank\tlineage\tuniversal_preference\tbest_ssu_cm\tcm_concordance\tassignment_status\tcluster_key\tassigned_model\tquantum_identity\ttrace_identity\tpaired_columns\tmember_paired_fraction\tmodel_paired_fraction");
			for(Member m:members){final int bucket=m.model<0?models.size():m.model;final String label=m.model<0?"UNASSIGNED":models.get(m.model).id;
				b.clear().append(arm).tab().append(m.id).tab().append(m.uniqueId).tab().append(m.taxid).tab().append(m.sequence.bases.length).tab().append(m.partial).tab().append(m.pseudo).tab().append(m.product).tab().append(binName(m.bin)).tab().append(rankName).tab().append(m.lineage).tab().append(m.preference).tab().append(m.cm).tab().append(m.concordance).tab().append(m.status).tab().append(m.clusterKey).tab().append(label);
				if(m.task==null){b.append("\tNA\tNA\tNA\tNA\tNA");}else{final NcrnaMitoSsuCompare.Result r=m.task.result;b.tab().append(r.identity,9).tab().append(r.traceIdentity,9).tab().append(r.paired).tab().append((double)r.paired/m.sequence.bases.length,9).tab().append((double)r.paired/models.get(m.model).length(),9);}
				Summary[] lineage=perLineage.get(m.lineage);if(lineage==null){lineage=summaries();perLineage.put(m.lineage,lineage);lineageOrder.add(m.lineage);}
				raw.print(b.nl());for(int s=0;s<STRATA.length;s++){if(inStratum(m,s)){overall[s].add(m);perModel[bucket][s].add(m);lineage[s].add(m);}}
				comp[bucket].add(m);
			}
			pairs.println("arm\tkind\tleft_id\tright_id\tleft_length\tright_length\tquery_side\tquantum_identity\ttrace_identity\tref_start0\tref_stop0\tpaired_columns\tleft_paired_fraction\tright_paired_fraction");
			for(Task t:tasks){if(t.kind.equals("MEMBER")){continue;}final NcrnaMitoSsuCompare.Result r=t.result;
				pairs.print(b.clear().append(arm).tab().append(t.kind).tab().append(t.left).tab().append(t.right).tab().append(t.a.length).tab().append(t.b.length).tab().append(r.candidateQuery?"left":"right").tab().append(r.identity,9).tab().append(r.traceIdentity,9).tab().append(r.start).tab().append(r.stop).tab().append(r.paired).tab().append((double)r.paired/t.a.length,9).tab().append((double)r.paired/t.b.length,9).nl());
			}
			summary.println("arm\tscope\tgroup\tstratum\tarm_copies\tcopies\tassigned_copies\tunassigned_copies\tunique_sequences\tunique_with_assignment\tunique_with_unassigned_copy\tcopy_mean_identity\tunique_mean_identity\tunique_weighting\tlineage_rank");
			for(int s=0;s<STRATA.length;s++){writeSummary(summary,b,"ARM","ALL",s,overall[s]);}
			for(int i=0;i<perModel.length;i++){final String label=i<models.size()?models.get(i).id:"UNASSIGNED";for(int s=0;s<STRATA.length;s++){writeSummary(summary,b,"CLUSTER",label,s,perModel[i][s]);}}
			int lineageCopies=0;for(String label:lineageOrder){final Summary[] rows=perLineage.get(label);lineageCopies+=rows[0].copies;for(int s=0;s<STRATA.length;s++){writeSummary(summary,b,"LINEAGE",label,s,rows[s]);}}
			require(lineageCopies==expected,"Lineage rows, including taxonomy gaps, must partition the complete arm");
			composition.println("arm\tmodel\tdimension\tvalue\tlength_bin\tcopies\tunique_sequences\tarm_copies");
			for(int i=0;i<comp.length;i++){final String label=i<models.size()?models.get(i).id:"UNASSIGNED";for(Cell c:comp[i].cells){composition.print(b.clear().append(arm).tab().append(label).tab().append(c.dimension).tab().append(c.value).tab().append(c.bin<0?"NA":binName(c.bin)).tab().append(c.copies).tab().append(c.seen.size()).tab().append(expected).nl());}}
			require(overall[0].copies==expected && overall[0].assigned==assigned,"Summary must conserve every original member including unassigned");
		}finally{RrnaResourceIO.finish(raw,pairs,summary,composition);}
	}
	private void writeSummary(ByteStreamWriter out,ByteBuilder b,String scope,String label,int stratum,Summary s){
		assert(s.assigned<=s.copies):"Summary denominators separate assigned from all copies";double uniqueSum=0;int uniqueAssigned=0,uniqueUnassigned=0;
		for(UniqueMean u:s.values){if(u.assigned>0){uniqueSum+=u.sum/u.assigned;uniqueAssigned++;}if(u.unassigned>0){uniqueUnassigned++;}}
		b.clear().append(arm).tab().append(scope).tab().append(label).tab().append(STRATA[stratum]).tab().append(expected).tab().append(s.copies).tab().append(s.assigned).tab().append(s.copies-s.assigned).tab().append(s.values.size()).tab().append(uniqueAssigned).tab().append(uniqueUnassigned).tab();
		if(s.assigned==0){b.append("NA");}else{b.append(s.sum/s.assigned,9);}b.tab();if(uniqueAssigned==0){b.append("NA");}else{b.append(uniqueSum/uniqueAssigned,9);}
		out.print(b.append("\tMEAN_OF_WITHIN_SEQUENCE_ASSIGNED_COPY_MEANS").tab().append(scope.equals("LINEAGE")?rankName:"NA").nl());
	}
	private String binName(int i){assert(i>=0 && i<=bins.length):"Length bin comes from validated ascending boundaries";return (i==0?1:bins[i-1])+"-"+(i==bins.length?"INF":Integer.toString(bins[i]-1));}
	private static boolean inStratum(Member m,int s){return s==0 || s==1 && !m.partial && !m.pseudo || s==2 && m.partial || s==3 && m.pseudo;}
	private static Summary[] summaries(){final Summary[] s=new Summary[STRATA.length];for(int i=0;i<s.length;i++){s[i]=new Summary();}return s;}
	private static boolean flag(String s){require(s.equals("true") || s.equals("false"),"Source flags must be literal booleans");return s.equals("true");}
	private static void require(boolean b,String message){if(!b){throw new IllegalArgumentException(message);}}
	@Override public void accumulate(Worker w){success&=w.ok;if(w.failure!=null){w.failure.printStackTrace(System.err);}}
	@Override public boolean success(){return success;}
	@Override public ReadWriteLock rwlock(){return null;}
	static final class Worker extends Thread{
		Worker(NcrnaMitoClusterEvaluation x){owner=x;setDaemon(true);}
		@Override public void run(){try{final AlignmentStats stats=new AlignmentStats(true);for(int start;(start=owner.next.getAndAdd(8))<owner.tasks.size();){final int end=Math.min(start+8,owner.tasks.size());for(int i=start;i<end;i++){final Task t=owner.tasks.get(i);t.result=NcrnaMitoSsuCompare.compare(t.a,t.b,stats);}}ok=true;}catch(Throwable t){failure=t;}}
		final NcrnaMitoClusterEvaluation owner;boolean ok;Throwable failure;
	}
	private static final class Task{Task(String k,String l,String r,byte[] aa,byte[] bb){kind=k;left=l;right=r;a=aa;b=bb;}final String kind,left,right;final byte[] a,b;NcrnaMitoSsuCompare.Result result;}
	private static final class Sequence{Sequence(String id,byte[] b,int i){uniqueId=id;bases=b;index=i;}final String uniqueId;final byte[] bases;final int index;final ObjectMap<String,Task> assigned=new ObjectMap<String,Task>(String.class,Task.class);}
	private static final class Member{String id,header,uniqueId,product,preference,cm,concordance,status,clusterKey,lineage;int taxid,bin,model=-1;boolean partial,pseudo,assignmentSeen,auditSeen;Sequence sequence;Task task;}
	private static final class UniqueMean{int assigned,unassigned;double sum;}
	private static final class Summary{
		void add(Member m){copies++;int index=ids.get(m.sequence.index);if(index<0){index=values.size();ids.put(m.sequence.index,index);values.add(new UniqueMean());}final UniqueMean u=values.get(index);if(m.task==null){u.unassigned++;}else{final double v=m.task.result.identity;assigned++;sum+=v;u.assigned++;u.sum+=v;}}
		int copies,assigned;double sum;final IntHashMap ids=new IntHashMap();final ArrayList<UniqueMean> values=new ArrayList<UniqueMean>();
	}
	private static final class Cell{Cell(String d,String v,int b){dimension=d;value=v;bin=b;}void add(Member m){copies++;seen.put(m.sequence.index,1);}final String dimension,value;final int bin;int copies;final IntHashMap seen=new IntHashMap();}
	private static final class Composition{
		Composition(int n){binCount=n;}
		void add(Member m){Cell[] bins=products.get(m.product);if(bins==null){bins=new Cell[binCount];products.put(m.product,bins);}Cell c=bins[m.bin];if(c==null){c=new Cell("PRODUCT_LENGTH",m.product,m.bin);bins[m.bin]=c;cells.add(c);}c.add(m);
			add(lineages,"LINEAGE",m.lineage,m);add(preferences,"UNIVERSAL_PREFERENCE",m.preference,m);add(cms,"BEST_SSU_CM",m.cm,m);add(statuses,"ASSIGNMENT_STATUS",m.status,m);
		}
		void add(ObjectMap<String,Cell> map,String dim,String value,Member m){Cell c=map.get(value);if(c==null){c=new Cell(dim,value,-1);map.put(value,c);cells.add(c);}c.add(m);}
		final int binCount;final ObjectMap<String,Cell[]> products=new ObjectMap<String,Cell[]>(String.class,Cell[].class);
		final ObjectMap<String,Cell> lineages=new ObjectMap<String,Cell>(String.class,Cell.class),preferences=new ObjectMap<String,Cell>(String.class,Cell.class),cms=new ObjectMap<String,Cell>(String.class,Cell.class),statuses=new ObjectMap<String,Cell>(String.class,Cell.class);final ArrayList<Cell> cells=new ArrayList<Cell>();
	}
	private String input,assignments,consensus,hbm,audit,treeFile,rankName,arm,prefix;private int expected,rank,assigned;private boolean success=true;
	private int[] bins={200,500,1000,1500,2000,3000,5000};private final String[] refFiles=new String[4],refIds=new String[4];private final byte[][] refs=new byte[4][];
	private final ObjectMap<String,Member> byId=new ObjectMap<String,Member>(String.class,Member.class);private final ObjectMap<String,Sequence> byUnique=new ObjectMap<String,Sequence>(String.class,Sequence.class);
	private final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class),modelTokens=new ObjectIntMap<String>(String.class);private final ArrayList<Member> members=new ArrayList<Member>();private final ArrayList<Sequence> sequences=new ArrayList<Sequence>();private final ArrayList<Read> models=new ArrayList<Read>();private final ArrayList<Task> tasks=new ArrayList<Task>();private final AtomicInteger next=new AtomicInteger();
	private static final String[] STRATA={"ALL","SOURCE_NONPARTIAL_NONPSEUDO","PARTIAL","PSEUDO"},REF_NAMES={"16S","18S","23S","eukLSU"};
	private static final String[] SUFFIXES={".members.tsv.gz",".pairs.tsv.gz",".summary.tsv",".composition.tsv"};
}
