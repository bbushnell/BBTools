package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMap;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

/** Exhaustive best-library identity, retaining every input copy and every score.
 * Uses the identical shorter-query Quantum/trace path as the cluster evaluator.
 * No seed, shortlist, identity exclusion or early acceptance is applied.
 * @author Brian, Raiden */
public final class NcrnaMitoCrossLibrary implements Accumulator<NcrnaMitoCrossLibrary.Worker> {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoCrossLibrary.class,false);
		try{final NcrnaMitoCrossLibrary x=new NcrnaMitoCrossLibrary();x.parse(pp.args);x.process();
			pp.outstream.println("MITO_CROSS_LIBRARY_PASS arm="+x.arm+" copies="+x.members.size()+" unique="+x.sequences.size()+" models="+x.models.size()+" alignments="+((long)x.sequences.size()*x.models.size())+" exhaustive=true no_exclusion=true");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();String verified=null;
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("ref")){reference=v;}else if(k.equals("tree")){treeFile=v;}else if(k.equals("outprefix")){prefix=v;}else if(k.equals("arm")){arm=v;}
			else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}else if(k.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(v);}else if(k.equals("verified")){verified=v;}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		input=p.in1;require(input!=null && reference!=null && treeFile!=null && prefix!=null && arm!=null && expected>0 && expectedModels>0 && !p.append && p.maxReads<0,"Require exact input/model denominators and fresh output; no append/read limit");
		require(verified!=null && treeFile.indexOf(',')<0,"Use explicitly pinned inputs and one taxonomy tree");final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);
		for(String file:new String[]{input,reference,treeFile}){require(pins.get(file)!=null,"Missing producer-verified input: "+file);}
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));RrnaResourceIO.fresh(prefix,SUFFIXES);
	}
	private void process(){
		final ObjectIntMap<String> loci=new ObjectIntMap<String>(String.class);final ObjectMap<String,Sequence> cache=new ObjectMap<String,Sequence>(String.class,Sequence.class);
		RrnaResourceIO.stream(input,r->{final String id=RrnaResourceIO.first(r.id),uid=RrnaResourceIO.field(r.id,"unique");require(loci.get(id)<0 && r.bases!=null && r.length()>0,"Duplicate locus or empty sequence: "+id);loci.put(id,1);
			Sequence s=cache.get(uid);if(s==null){s=new Sequence(uid,r.bases.clone(),sequences.size());cache.put(uid,s);sequences.add(s);}else{require(Arrays.equals(s.bases,r.bases),"Unique ID binds conflicting bases: "+uid);}
			final int tid=TaxTree.parseTaxID(r.id);require(tid>0,"Each source copy needs its actual host taxid");members.add(new Member(id,s,tid,flag(RrnaResourceIO.field(r.id,"partial")),flag(RrnaResourceIO.field(r.id,"pseudo"))));
		});require(members.size()==expected,"Full input denominator differs from expected copies");
		final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class);RrnaResourceIO.stream(reference,r->{final String token=RrnaResourceIO.first(r.id);require(modelIds.get(token)<0 && r.bases!=null && r.length()>0,"Empty or duplicate model token");modelIds.put(token,models.size());models.add(new Read(r.bases.clone(),null,r.id,models.size()));});
		require(models.size()==expectedModels,"The complete declared model library is required");
		final TaxTree tree;final boolean old=TaxTree.LOAD_FAKE_TREE;try{TaxTree.LOAD_FAKE_TREE=false;tree=TaxTree.loadTaxTreeFile(treeFile);}finally{TaxTree.LOAD_FAKE_TREE=old;}
		final IntHashMap tids=new IntHashMap();final ArrayList<String> lineages=new ArrayList<String>();
		for(Member m:members){int i=tids.get(m.taxid);if(i<0){i=lineages.size();lineages.add(NcrnaMitoClusterEvaluation.lineage(tree,m.taxid,TaxTree.PHYLUM_E));tids.put(m.taxid,i);}m.phylum=lineages.get(i);}
		final ArrayList<Worker> workers=new ArrayList<Worker>();for(int i=0;i<Math.min(Shared.threads(),sequences.size());i++){workers.add(new Worker(this));}
		require(ThreadWaiter.startAndWait(workers,this),"An exhaustive-alignment worker failed; incomplete model searches cannot become coverage evidence");write();
	}
	private void write(){
		assert(!sequences.isEmpty() && !models.isEmpty()):"Output follows complete input loading and successful exhaustive workers";
		final ByteStreamWriter raw=RrnaResourceIO.writer(prefix+SUFFIXES[0]),scores=RrnaResourceIO.writer(prefix+SUFFIXES[1]),summary=RrnaResourceIO.writer(prefix+SUFFIXES[2]),hist=RrnaResourceIO.writer(prefix+SUFFIXES[3]);
		final ByteBuilder b=new ByteBuilder();final Summary all=new Summary();final ObjectMap<String,Summary> groups=new ObjectMap<String,Summary>(String.class,Summary.class);final ArrayList<String> order=new ArrayList<String>();
		try{
			raw.println("arm\tlocus_id\tunique_id\ttaxid\tphylum\tpartial\tpseudo\tlength\tbest_model\tbest_identity\ttied_models\tmember_paired_fraction\tmodel_paired_fraction");
			for(Member m:members){final Sequence s=m.sequence;final Read model=models.get(s.best);final NcrnaMitoSsuCompare.Result r=s.result;require(r!=null && s.scores.length==models.size(),"Every sequence needs the complete model-score array and a best alignment");
				raw.print(b.clear().append(arm).tab().append(m.id).tab().append(s.id).tab().append(m.taxid).tab().append(m.phylum).tab().append(m.partial).tab().append(m.pseudo).tab().append(s.bases.length).tab().append(model.id).tab().append(r.identity,9).tab().append(s.ties).tab().append((double)r.paired/s.bases.length,9).tab().append((double)r.paired/model.length(),9).nl());
				Summary g=groups.get(m.phylum);if(g==null){g=new Summary();groups.put(m.phylum,g);order.add(m.phylum);}all.add(m);g.add(m);
			}
			scores.println("arm\tunique_id\tmodel\tquantum_identity");for(Sequence s:sequences){for(int i=0;i<models.size();i++){scores.print(b.clear().append(arm).tab().append(s.id).tab().append(models.get(i).id).tab().append(s.scores[i],9).nl());}}
			summary.println("arm\tscope\tphylum\tcopies\tunique_sequences\tpartial_copies\tpseudo_copies\tcopy_mean\tunique_mean\tcopy_ge_070\tcopy_ge_075\tcopy_ge_080\tunique_ge_070\tunique_ge_075\tunique_ge_080\tminimum\tmaximum\tmodels\tcomparison");
			hist.println("arm\tscope\tphylum\tbin_low\tbin_high_exclusive\tcopies\tunique_sequences\tbin_rule");
			writeGroup(summary,hist,b,"ALL","ALL",all);int copies=0;for(String name:order){final Summary g=groups.get(name);copies+=g.copies;writeGroup(summary,hist,b,"PHYLUM",name,g);}require(copies==expected && all.copies==expected && all.seen.size()==sequences.size(),"Phylum rows, including taxonomy gaps, must conserve the full copy denominator");
		}finally{RrnaResourceIO.finish(raw,scores,summary,hist);}
	}
	private void writeGroup(ByteStreamWriter out,ByteStreamWriter hist,ByteBuilder b,String scope,String name,Summary g){
		assert(g.copies>0 && g.seen.size()>0):"Every emitted summary represents retained source members";
		b.clear().append(arm).tab().append(scope).tab().append(name).tab().append(g.copies).tab().append(g.seen.size()).tab().append(g.partial).tab().append(g.pseudo).tab().append(g.sum/g.copies,9).tab().append(g.uniqueSum/g.seen.size(),9);
		for(int n:g.copyThresholds){b.tab().append(n);}for(int n:g.uniqueThresholds){b.tab().append(n);}out.print(b.tab().append(g.min,9).tab().append(g.max,9).tab().append(models.size()).append("\tMAXIMUM_OVER_ALL_MODELS_FIRST_FASTA_TIE").nl());
		int count=0,unique=0;for(int i=0;i<=100;i++){count+=g.hist[i];unique+=g.uniqueHist[i];b.clear().append(arm).tab().append(scope).tab().append(name).tab().append(i/100.0,2).tab();if(i==100){b.append("INCLUSIVE_1");}else{b.append((i+1)/100.0,2);}hist.print(b.tab().append(g.hist[i]).tab().append(g.uniqueHist[i]).append("\tFLOOR_NATIVE_FLOAT_IDENTITY_TIMES_100").nl());}
		require(count==g.copies && unique==g.seen.size(),"Identity histogram must conserve the matching summary denominators");
	}
	@Override public void accumulate(Worker w){ok&=w.ok;if(w.failure!=null){w.failure.printStackTrace(System.err);}}
	@Override public boolean success(){return ok;}
	@Override public ReadWriteLock rwlock(){return null;}
	static final class Worker extends Thread{
		Worker(NcrnaMitoCrossLibrary x){owner=x;setDaemon(true);}
		@Override public void run(){try{final AlignmentStats stats=new AlignmentStats(true);for(int start;!owner.abort && (start=owner.next.getAndAdd(4))<owner.sequences.size();){for(int i=start;i<Math.min(start+4,owner.sequences.size());i++){
			final Sequence s=owner.sequences.get(i);s.scores=new float[owner.models.size()];
			for(int m=0;m<owner.models.size();m++){final NcrnaMitoSsuCompare.Result r=NcrnaMitoSsuCompare.compare(s.bases,owner.models.get(m).bases,stats);s.scores[m]=r.identity;if(s.result==null || r.identity>s.result.identity){s.best=m;s.result=r;s.ties=1;}else if(r.identity==s.result.identity){s.ties++;}}
		}}ok=!owner.abort;}catch(Throwable t){failure=t;owner.abort=true;}}
		final NcrnaMitoCrossLibrary owner;boolean ok;Throwable failure;
	}
	private static final class Sequence{Sequence(String key,byte[] b,int i){id=key;bases=b;index=i;}final String id;final byte[] bases;final int index;float[] scores;int best=-1,ties;NcrnaMitoSsuCompare.Result result;}
	private static final class Member{Member(String key,Sequence s,int tid,boolean p,boolean ps){id=key;sequence=s;taxid=tid;partial=p;pseudo=ps;}final String id;final Sequence sequence;final int taxid;final boolean partial,pseudo;String phylum;}
	private static final class Summary{
		void add(Member m){final float value=m.sequence.result.identity;require(Float.isFinite(value) && value>=0 && value<=1,"Only a complete finite best alignment may enter a coverage summary");copies++;sum+=value;partial+=m.partial?1:0;pseudo+=m.pseudo?1:0;min=Math.min(min,value);max=Math.max(max,value);final int bin=Math.min(100,(int)(value*100));hist[bin]++;
			final boolean first=seen.get(m.sequence.index)<0;if(first){seen.put(m.sequence.index,1);uniqueSum+=value;uniqueHist[bin]++;}final int bits=thresholdBits(value);for(int i=0;i<THRESHOLDS.length;i++){if((bits&(1<<i))!=0){copyThresholds[i]++;if(first){uniqueThresholds[i]++;}}}}
		int copies,partial,pseudo;double sum,uniqueSum;float min=1,max=0;final int[] copyThresholds=new int[3],uniqueThresholds=new int[3],hist=new int[101],uniqueHist=new int[101];final IntHashMap seen=new IntHashMap();
	}
	private static boolean flag(String value){require(value.equals("true") || value.equals("false"),"Partial/pseudo source flags must be literal booleans");return value.equals("true");}
	/** Inclusive thresholds in the native float identity domain, lowest cut in bit0. */
	public static int thresholdBits(float value){require(Float.isFinite(value) && value>=0 && value<=1,"Threshold counts require a finite native identity in[0,1]");int bits=0;for(int i=0;i<THRESHOLDS.length;i++){if(value>=THRESHOLDS[i]){bits|=1<<i;}}return bits;}
	private static void require(boolean value,String message){if(!value){throw new IllegalArgumentException(message);}}
	private String input,reference,treeFile,prefix,arm;private int expected,expectedModels;private boolean ok=true;private volatile boolean abort;
	private final ArrayList<Member> members=new ArrayList<Member>();private final ArrayList<Sequence> sequences=new ArrayList<Sequence>();private final ArrayList<Read> models=new ArrayList<Read>();private final AtomicInteger next=new AtomicInteger();
	private static final float[] THRESHOLDS={.70f,.75f,.80f};private static final String[] SUFFIXES={".members.tsv.gz",".scores.tsv.gz",".summary.tsv",".histogram.tsv"};
}
