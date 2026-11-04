package test.prok;

import prok.MitoSsuModelPolicy;
import prok.RrnaResourceIO;
import prok.TrnaConsensusBuilder;
import prok.TrnaKmerIndex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import tax.TaxTree;
import template.ThreadWaiter;

/** Construction assignments, not a caller or an all-pairs coverage matrix.
 * Rank the153 frozen mitochondrial models by shared k9 with each exact source
 * gene, align at most8, retain the maximum shorter-query identity. Then align
 * that one model to the original genomic flank for raw endpoint characterization.
 * Every input member and its original pilot caller status remain in the ledger.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointAssignments {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,MitoEndpointAssignments.class,false);
		try{final MitoEndpointAssignments x=new MitoEndpointAssignments();x.parse(pp.args);x.run();}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("ledger")){ledger=v;}else if(k.equals("ssu")){ssu=v;}else if(k.equals("insect")){insect=v;}
			else if(k.equals("pilotanchors")){pilotAnchors=v;}else if(k.equals("consensus")){consensus=v;}
			else if(k.equals("outprefix")){prefix=v;}else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}
			else if(k.equals("passes")){passes=Integer.parseInt(v);}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}}
		fasta=p.in1;require(fasta!=null && ledger!=null && ssu!=null && insect!=null && consensus!=null && pilotAnchors!=null && prefix!=null,
			"Assignments need the complete source/window roster, both exact-source arms, shipped library and original pilot statuses");
		require((expected==256 || expected==16110) && !p.append && p.maxReads<0,"Only the declared pilot or whole construction population may be processed");
		require(passes>=1 && passes<=8 && (passes==1 || expected==256 && Shared.threads()==1),"Repeated timing is confined to the256member single-thread pilot; full assignment runs once");
		Shared.TRIM_READ_DESCRIPTION=false;stream.Read.TO_UPPER_CASE=true;
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
		RrnaResourceIO.fresh(prefix,".members.tsv.gz",".summary.tsv.gz");
	}
	private void run(){
		final byte[][] all=TrnaConsensusBuilder.loadLibrary(consensus);names=TrnaConsensusBuilder.lastLoadedNames;
		require(all.length==158 && MitoSsuModelPolicy.coverage("16S",names),"Construction mapping must use the unchanged153mitochondrial models after the original five");
		library=Arrays.copyOfRange(all,5,all.length);load();
		for(int pass=0;pass<passes;pass++){
			final AtomicInteger next=new AtomicInteger();final ArrayList<Worker> workers=new ArrayList<Worker>();
			for(int i=0;i<Math.min(Shared.threads(),rows.size());i++){workers.add(new Worker(rows,library,next));}
			require(ThreadWaiter.startAndWait(workers),"All construction-mapping workers must terminate before publication");
			for(Worker w:workers){if(w.failure!=null){throw new IllegalStateException("Construction mapping failed",w.failure);}require(w.success,"Incomplete construction worker");}
			for(Row r:rows){
				if(pass==0){r.baseline=new int[]{r.selection.model,Float.floatToIntBits(r.selection.identity),r.selection.tried,r.selection.hits,r.start,r.stop,Float.floatToIntBits(r.identity),r.state};}
				else{require(Arrays.equals(r.baseline,new int[]{r.selection.model,Float.floatToIntBits(r.selection.identity),r.selection.tried,r.selection.hits,r.start,r.stop,Float.floatToIntBits(r.identity),r.state}),
					"Timing repetitions must reproduce every member assignment, endpoint and disposition exactly");}
			}
		}
		write();
	}
	private void load(){
		final ByteFile in=RrnaResourceIO.open(ledger);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,MitoEndpointCorpus.HEADER,ledger);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18 && p.termEquals("true",12),"Assignments consume a complete selected-source roster");
				final Row r=new Row();r.id=p.parseString(0);r.arm=p.parseString(1);r.host=p.parseString(3);r.split=p.parseString(4);r.exactLength=p.parseInt(5);
				r.partial=flag(p,9);r.pseudo=flag(p,10);r.source=p.parseString(13);r.ts=p.parseInt(14);r.te=p.parseInt(15);r.windowLength=p.parseInt(16);
				require((r.arm.equals("SSU") || r.arm.equals("INSECT")) && (r.split.equals("TRAIN") || r.split.equals("VALIDATION")) && r.exactLength>0 && map.put(r.id,r)==null,
					"Every source must retain a unique ID, source arm, host split and positive length");rows.add(r);
				final String old=hostSplits.put(r.host,r.split);require(old==null || old.equals(r.split),"One host cannot enter both development splits");
			}
		}finally{require(!in.close(),"Source roster read failed");}
		require(rows.size()==expected,"Construction mapping must conserve the full declared member population");
		loadExact(ssu,"SSU",11420);loadExact(insect,"INSECT",7050);
		RrnaResourceIO.stream(fasta,read->{final Row r=map.get(RrnaResourceIO.first(read.id));require(r!=null && r.window==null && read.length()==r.windowLength,"Each genomic flank must bind one source member");
			require(Integer.parseInt(RrnaResourceIO.field(read.id,"start0"))==r.ts && Integer.parseInt(RrnaResourceIO.field(read.id,"stop0"))==r.te
				&& RrnaResourceIO.field(read.id,"host").equals(r.host) && RrnaResourceIO.field(read.id,"split").equals(r.split),"Flank coordinates and host split must remain the frozen source values");r.window=read.bases;});
		for(Row r:rows){require(r.exact!=null,"Every declared construction member needs its original exact gene: "+r.id);
			require(r.source.equals("UNSUPPORTED_WINDOW")?r.window==null:r.window!=null && r.ts>=0 && r.te>=r.ts && r.te<r.window.length,"Explicit source geometry must account for every member");}
		final ByteFile pilot=RrnaResourceIO.open(pilotAnchors);int n=0;
		try{final byte[] header=pilot.nextLine();require(header!=null,"Missing pilot caller ledger");p.set(header);
			require(p.terms()==18 && p.termEquals("locus_id",0) && p.termEquals("anchor_status",5),"Pilot caller schema changed");
			for(byte[] line;(line=pilot.nextLine())!=null;){p.set(line);require(p.terms()==18,"Pilot row width changed");final Row r=map.get(p.parseString(0));
				require(r!=null && r.pilotStatus.equals("NOT_IN_PILOT") && p.termEquals(r.host,2) && p.termEquals(r.split,3),"All256original pilot rows must survive construction enrollment exactly once");
				r.pilotStatus=p.parseString(5);n++;}
		}finally{require(!pilot.close(),"Pilot status ledger read failed");}
		require(n==256,"Keep all original256pilot members, including the58caller-noneligible members");
	}
	private void loadExact(String path,String arm,int expectedInput){
		final long[] count={0};RrnaResourceIO.stream(path,read->{count[0]++;final Row r=map.get(RrnaResourceIO.first(read.id));if(r==null){return;}
			require(r.exact==null && r.arm.equals(arm) && read.length()==r.exactLength && TaxTree.parseTaxID(read.id)==Integer.parseInt(r.host)
				&& RrnaResourceIO.field(read.id,"partial").equals(Boolean.toString(r.partial)) && RrnaResourceIO.field(read.id,"pseudo").equals(Boolean.toString(r.pseudo)),
				"Exact gene must bind the sealed member length, host, arm and source flags");r.exact=read.bases;});
		require(count[0]==expectedInput,"Use the complete archived source arm, not a recall-selected subset");
	}
	/** The shortlist is bounded and no identity threshold silently excludes a gene. */
	public static Selection assign(byte[] exact,byte[][] library,TrnaKmerIndex index){
		require(exact!=null && exact.length>0 && library.length>0 && index!=null,"Construction mapping requires real gene/library bases and the shared index");
		final int[] shortlist=index.shortlist(exact,8),pos=new int[2];final Selection best=new Selection();best.tried=shortlist.length;
		require(shortlist.length<=8,"G11's bounded construction mapping never expands to an all-pairs matrix");
		for(int m:shortlist){final byte[] model=library[m];
			final float id=exact.length<=model.length?alignIndependent(exact,model,pos):alignIndependent(model,exact,pos);
			require(Float.isFinite(id) && id>=0 && id<=1,"Assignment identities must be finite shorter-query Quantum scores");
			if(best.model<0 || id>best.identity){best.model=m;best.identity=id;best.hits=index.lastSharedCount(m);}
		}
		return best;
	}
	/** pos[0] is also an input hint to IDAlignerStatics.decideBandwidth. Reset it
	 * to the fresh-array baseline so model/locus scheduling cannot alter the fit. */
	public static float alignIndependent(byte[] query,byte[] reference,int[] pos){
		require(query!=null && reference!=null && query.length>0 && query.length<=reference.length && pos!=null && pos.length>=2,
			"Independent full-query fits require valid sequence lengths and reference-coordinate storage");
		pos[0]=pos[1]=0;return QuantumAligner.alignStatic(query,reference,pos);
	}
	private void write(){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".members.tsv.gz"),summary=RrnaResourceIO.writer(prefix+".summary.tsv.gz");final ByteBuilder b=new ByteBuilder();
		final long[][] sums=new long[STATES.length][4];long mappings=0,genomic=0;
		try{out.println("locus_id\tarm\thost\tsplit\tsource_status\tpilot_caller_status\tdisposition\tmodel_index\tmodel\tshared9\tmodels_aligned\tassignment_identity_shorter_query\tgenomic_identity_full_model\ttruth_start0\ttruth_stop0\traw_start0\traw_stop0\tcorrection5\tcorrection3\tpartial\tpseudo");
			for(Row r:rows){require(r.selection!=null,"Every original member requires a terminal assignment disposition");final Selection s=r.selection;final boolean fitted=r.start>=0;
				b.clear().append(r.id).tab().append(r.arm).tab().append(r.host).tab().append(r.split).tab().append(r.source).tab().append(r.pilotStatus).tab().append(STATES[r.state])
					.tab().append(s.model<0?-1:s.model+5).tab().append(s.model<0?"NA":names[s.model+5]).tab().append(s.hits).tab().append(s.tried).tab();
				if(s.model<0){b.append("NA");}else{b.appendSlow(s.identity);}b.tab();if(!fitted){b.append("NA");}else{b.appendSlow(r.identity);}
				b.tab().append(r.ts).tab().append(r.te).tab().append(r.start).tab().append(r.stop).tab();
				if(fitted){b.append(r.ts-r.start).tab().append(r.te-r.stop);}else{b.append("NA\tNA");}
				out.print(b.tab().append(r.partial).tab().append(r.pseudo).nl());
				final long[] c=sums[r.state];c[0]++;c[r.split.equals("TRAIN")?1:2]++;if(!r.pilotStatus.equals("NOT_IN_PILOT")){c[3]++;}mappings+=s.tried;genomic+=fitted?1:0;
			}
			summary.println("disposition\tmembers\ttrain\tvalidation\tpilot_members\tassignment_alignments\tgenomic_alignments\tclaim");long total=0,pilot=0;
			for(int i=0;i<STATES.length;i++){final long[] c=sums[i];total+=c[0];pilot+=c[3];require(c[0]==c[1]+c[2],"Disposition splits must conserve every member");
				summary.print(b.clear().append(STATES[i]).tab().append(c[0]).tab().append(c[1]).tab().append(c[2]).tab().append(c[3]).append("\tNA\tNA\tCONSTRUCTION_NOT_CALLER_ACCURACY").nl());}
			require(total==expected && pilot==256,"Assignment rows and dispositions must preserve construction and original pilot denominators");
			summary.print(b.clear().append("ALL").tab().append(total).append("\tNA\tNA\t").append(pilot).tab().append(mappings).tab().append(genomic).append("\tTOP8_K9_MAX_IDENTITY_NOT_ALL_PAIRS").nl());
		}finally{RrnaResourceIO.finish(out,summary);}
		System.err.println("MITO_ENDPOINT_ASSIGNMENTS_PASS members="+rows.size()+" original_pilot_members=256 assignment_alignments="+mappings+" genomic_alignments="+genomic+" timing_passes="+passes+" actual_alignment_calls="+((mappings+genomic)*passes)+" models=153 topn=8 k=9 no_fit=true");
	}
	private static final class Worker extends Thread{
		Worker(ArrayList<Row> r,byte[][] l,AtomicInteger n){rows=r;library=l;next=n;index=new TrnaKmerIndex(l,9,false,1,0,0,1);}
		@Override public void run(){try{final int[] pos=new int[2];for(int batch;(batch=next.getAndAdd(4))<rows.size();){for(int i=batch;i<Math.min(batch+4,rows.size());i++){
			final Row r=rows.get(i);final Selection s=r.selection=assign(r.exact,library,index);
			if(s.model<0){r.state=0;continue;}if(r.window==null){r.state=1;continue;}
			if(library[s.model].length>r.window.length){r.state=2;continue;}
			r.identity=alignIndependent(library[s.model],r.window,pos);r.start=pos[0];r.stop=pos[1];
			require(Float.isFinite(r.identity) && r.identity>=0 && r.identity<=1 && r.start>=0 && r.stop>=r.start && r.stop<r.window.length,"Genomic full-model fit must retain valid identity and inclusive reference endpoints");
			final int overlap=Math.max(0,Math.min(r.stop,r.te)-Math.max(r.start,r.ts)+1);
			r.state=r.partial||r.pseudo?3:2L*overlap<=r.te-r.ts+1L || 2L*overlap<=r.stop-r.start+1L?4:5;
		}}success=true;}catch(Throwable t){failure=t;next.set(rows.size());}}
		final ArrayList<Row> rows;final byte[][] library;final AtomicInteger next;final TrnaKmerIndex index;boolean success;Throwable failure;
	}
	public static final class Selection{public int model=-1,hits,tried;public float identity=Float.NaN;}
	private static final class Row{String id,arm,host,split,source,pilotStatus="NOT_IN_PILOT";byte[] exact,window;boolean partial,pseudo;int exactLength,windowLength,ts,te,start=-1,stop=-1,state;int[] baseline;float identity=Float.NaN;Selection selection;}
	private static boolean flag(LineParser1 p,int i){require(p.termEquals("true",i) || p.termEquals("false",i),"Source flags must remain explicit");return p.termEquals("true",i);}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final ArrayList<Row> rows=new ArrayList<Row>();private final ObjectMap<String,Row> map=new ObjectMap<String,Row>(String.class,Row.class);
	private final ObjectMap<String,String> hostSplits=new ObjectMap<String,String>(String.class,String.class);
	private String fasta,ledger,ssu,insect,consensus,pilotAnchors,prefix;private int expected,passes=1;private byte[][] library;private String[] names;
	private static final String[] STATES={"NO_SHARED_9MER_MODEL","UNSUPPORTED_WINDOW","MODEL_LONGER_THAN_WINDOW","UNCERTAIN_SOURCE_ENDS","GENOMIC_FIT_NOT_ASSOCIATED","COMPLETE_SOURCE_ASSOCIATED"};
}
