package test.prok;

import prok.CallGenes;
import prok.ConservedRnaSeedIndex;
import prok.GeneCaller;
import prok.MitoSsuTierPolicy;
import prok.NcrnaBoundaryInstrumentSink;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.Prok16sEndpointConfig;
import prok.RrnaResourceIO;

import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.LongHashSet;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Capture raw mitochondrial anchors using the shipped two-tier scavenger.
 * Input windows are source-oriented construction loci, not intact-genome accuracy
 * benchmarks. Misses, original-tier hits and ambiguous associations stay explicit.
 * Captured mitochondrial coordinates are raw; original-tier coordinates may
 * already have passed their shipped endpoint nets and are excluded from ranges.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointAnchors {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,MitoEndpointAnchors.class,false);
		try{final MitoEndpointAnchors x=new MitoEndpointAnchors();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("ledger")){ledger=v;}else if(k.equals("outprefix")){prefix=v;}else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}}
		fasta=p.in1;require(ledger!=null && fasta!=null && prefix!=null && expected>0 && !p.append && p.maxReads<0,"Complete pilot ledger/FASTA and exact expected population required");
		Shared.TRIM_READ_DESCRIPTION=false;stream.Read.TO_UPPER_CASE=true;ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
		RrnaResourceIO.fresh(prefix,".loci.tsv.gz",".calls.tsv.gz",".ranges.tsv.gz");
	}
	private void run(){
		load();CallGenes.setRrna17Profile("16s");CallGenes.loadRrna17Resources("16s");
		new Prok16sEndpointConfig().apply(GeneCaller.ncrnaFamilies);
		require(GeneCaller.ncrnaFamilies.size()==1,"Anchor capture must load exactly the shipped16S family");
		final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(0);final MitoSsuTierPolicy policy=MitoSsuTierPolicy.configured();
		require(family.library.length==158 && family.kLong==17 && family.kmerSet.size()==680 && family.indexK==9 && family.indexTopN==1
			&& family.seedMinHits==2 && !family.adaptive && family.fixedMinHits==1 && family.idPass==.66f && family.reuseConsensusAlignment
			&& !family.trimAlignmentExtent && family.preserveContainingWindows && family.modelThresholds!=null
			&& policy.maxModels==3 && policy.minHits==80 && policy.queryFraction==0 && policy.topFraction==0,
			"Pilot must preserve the shipped153+5 library,680 seeds and two-tier controls; overrides require a new comparison");
		final AtomicInteger next=new AtomicInteger();final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(Shared.threads(),rows.size()); i++){workers.add(new Worker(rows,next,family));}
		require(ThreadWaiter.startAndWait(workers),"Anchor workers must finish before publishing reachability");
		long alignments=0;for(Worker w:workers){if(w.failure!=null){throw new IllegalStateException("Anchor worker failed",w.failure);}require(w.success,"Incomplete anchor worker");alignments+=w.caller.alignmentCount();}
		write(family);System.err.println("MITO_ENDPOINT_ANCHORS_PASS members="+rows.size()+" alignments="+alignments+" output_files=3 scope=CONSTRUCTION_COST_PILOT no_training=true");
	}
	private void load(){
		final ObjectMap<String,Row> byId=new ObjectMap<String,Row>(String.class,Row.class);final ByteFile in=RrnaResourceIO.open(ledger);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,MitoEndpointCorpus.HEADER,ledger);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18 && p.termEquals("true",12),"Use the complete selected pilot ledger");
				final Row r=new Row();r.id=p.parseString(0);r.arm=p.parseString(1);r.host=p.parseString(3);r.split=p.parseString(4);r.source=p.parseString(13);r.start=p.parseInt(14);r.stop=p.parseInt(15);r.length=p.parseInt(16);r.gc=p.parseFloat(17);
				require(byId.put(r.id,r)==null && (r.split.equals("TRAIN") || r.split.equals("VALIDATION")),"Pilot source IDs/splits must remain unique and declared");rows.add(r);}
		}finally{require(!in.close(),"Pilot ledger read failed");}
		require(rows.size()==expected,"Pilot population changed");
		RrnaResourceIO.stream(fasta,read->{final Row r=byId.get(RrnaResourceIO.first(read.id));require(r!=null && r.bases==null && read.length()==r.length,"Each flank must bind one selected source row");
			require(Integer.parseInt(RrnaResourceIO.field(read.id,"start0"))==r.start && Integer.parseInt(RrnaResourceIO.field(read.id,"stop0"))==r.stop
				&& RrnaResourceIO.field(read.id,"host").equals(r.host) && RrnaResourceIO.field(read.id,"split").equals(r.split),"Flank metadata disagrees with frozen ledger");r.bases=read.bases;});
		for(Row r:rows){require(r.source.equals("UNSUPPORTED_WINDOW")?r.bases==null:r.bases!=null && r.start>=0 && r.stop>=r.start && r.stop<r.bases.length,"Explicit geometry status must account for every pilot member: "+r.id);}
	}
	private void write(NcrnaFamily f){
		final ByteStreamWriter loci=RrnaResourceIO.writer(prefix+".loci.tsv.gz"),calls=RrnaResourceIO.writer(prefix+".calls.tsv.gz"),ranges=RrnaResourceIO.writer(prefix+".ranges.tsv.gz");final ByteBuilder b=new ByteBuilder();
		try{
			loci.println("locus_id\tarm\thost\tsplit\tsource_status\tanchor_status\treturned_calls\toverlapping_calls\tmodel_index\tmodel\ttruth_start0\ttruth_stop0\tcaptured_start0\tcaptured_stop0\tcorrection5\tcorrection3\tidentity\tgc");
			calls.println("locus_id\tmodel_index\tmodel\twindow_start0\twindow_stop0\tcaptured_start0\tcaptured_stop0\tidentity\treciprocal_overlap");
			for(Row r:rows){
				Capture chosen=null;int overlap=0;
				for(Capture c:r.calls){final int n=Math.max(0,Math.min(c.stop,r.stop)-Math.max(c.start,r.start)+1);c.overlap=2L*n>r.stop-r.start+1L && 2L*n>c.stop-c.start+1L;
					if(c.overlap){overlap++;chosen=c;}
					calls.print(b.clear().append(r.id).tab().append(c.model).tab().append(f.modelNames[c.model]).tab().append(c.ws).tab().append(c.we).tab().append(c.start).tab().append(c.stop).tab().appendSlow(c.identity).tab().append(c.overlap).nl());}
				r.status=r.bases==null?"UNSUPPORTED_WINDOW":overlap==0?"NO_ASSOCIATED_CALL":overlap>1?"AMBIGUOUS_CALLS":chosen.model<5?"ORIGINAL_TIER":!r.source.equals("EXACT_SOURCE_ENDS")?"UNCERTAIN_ENDS":"MITO_EXACT_SOURCE_ANCHOR";
				r.chosen=overlap==1?chosen:null;
				b.clear().append(r.id).tab().append(r.arm).tab().append(r.host).tab().append(r.split).tab().append(r.source).tab().append(r.status).tab().append(r.calls.size()).tab().append(overlap);
				if(r.chosen==null){b.append("\t-1\tNA");}else{b.tab().append(chosen.model).tab().append(f.modelNames[chosen.model]);}
				b.tab().append(r.start).tab().append(r.stop);
				if(r.chosen==null){b.append("\tNA\tNA\tNA\tNA\tNA");}
				else{b.tab().append(chosen.start).tab().append(chosen.stop);if(r.status.equals("MITO_EXACT_SOURCE_ANCHOR")){b.tab().append(MitoEndpointGeometry.correction(r.start,chosen.start)).tab().append(MitoEndpointGeometry.correction(r.stop,chosen.stop));}else{b.append("\tNA\tNA");}b.tab().appendSlow(chosen.identity);}
				loci.print(b.tab().append(r.gc,8).nl());
			}
			ranges.println("split\tradius\tall_members\teligible_mito_anchors\tstart_reachable\tstop_reachable\tboth_reachable\tclaim");
			for(String split:new String[]{"TRAIN","VALIDATION","ALL"}){for(int radius:new int[]{0,4,8,16,32,64,128,256,512,1024}){
				int all=0,n=0,s=0,e=0,both=0;
				for(Row r:rows){if(!split.equals("ALL") && !split.equals(r.split)){continue;}all++;if(!r.status.equals("MITO_EXACT_SOURCE_ANCHOR")){continue;}n++;
					final boolean a=Math.abs(MitoEndpointGeometry.correction(r.start,r.chosen.start))<=radius,z=Math.abs(MitoEndpointGeometry.correction(r.stop,r.chosen.stop))<=radius;
					s+=a?1:0;e+=z?1:0;both+=a && z?1:0;}
				ranges.print(b.clear().append(split).tab().append(radius).tab().append(all).tab().append(n).tab().append(s).tab().append(e).tab().append(both).append("\tPILOT_REACHABILITY_NOT_NN_ACCURACY").nl());
			}}
		}finally{RrnaResourceIO.finish(loci,calls,ranges);}
	}
	private static NcrnaScavenger scavenger(NcrnaFamily f){
		assert(f.rrnaEndpointFeatures!=null):"Original models retain their shipped endpoint behavior";
		final NcrnaScavenger s=new NcrnaScavenger(f.library,f.models,f.modelNames,f.kmerSet,f.kLong,f.minLen,f.windowPad,
			f.indexK,f.indexTopN,f.adaptive,f.adaptFloor,f.adaptTopFrac,f.adaptQFrac,f.fixedMinHits,f.scoreA,f.scoreB,f.idPass,f.idBorderline,
			f.boundary5NetTemplate,f.boundary3NetTemplate,f.boundaryStartTable,f.boundaryStopTable,f.boundaryStartInside,f.boundaryStartOutside,
			f.boundaryStopInside,f.boundaryStopOutside,f.boundaryMeanLen,f.boundaryStartOffsets,f.boundaryStopOffsets,f.boundaryMarginStart,f.boundaryMarginStop);
		GeneCaller.applyNcrnaFamilyControls(s,f);return s;
	}
	private static final class Worker extends Thread implements NcrnaBoundaryInstrumentSink{
		Worker(ArrayList<Row> r,AtomicInteger n,NcrnaFamily f){rows=r;next=n;caller=scavenger(f);caller.setInstrumentSink(this);index=new ConservedRnaSeedIndex(new LongHashSet[]{f.kmerSet},new boolean[]{false});}
		@Override public void run(){try{for(int i;(i=next.getAndIncrement())<rows.size();){current=rows.get(i);if(current.bases==null){continue;}
			final ConservedRnaSeedIndex.ScanResult hits=index.scan(current.bases);
			final ArrayList<Orf> found=caller.scavenge(current.id,current.bases,0,new ArrayList<int[]>(),hits.hits(0),hits.keys(0));
			require(found.size()==current.calls.size(),"Captured anchors must conserve every returned scavenger call");}success=true;
		}catch(Throwable t){failure=t;next.set(rows.size());}}
		@Override public void capture(String contig,int strand,int model,int ws,int we,int start,int stop,byte[] copy,int offset,boolean trimmed,boolean nn,int source,float identity,int alignedLength){
			require(current!=null && current.id.equals(contig) && strand==0 && !trimmed && !nn && start>=0 && stop<current.bases.length && start<=stop && Float.isFinite(identity),"Capture must preserve actual accepted sense-oriented raw anchors");
			current.calls.add(new Capture(model,ws,we,start,stop,identity));
		}
		final ArrayList<Row> rows;final AtomicInteger next;final NcrnaScavenger caller;final ConservedRnaSeedIndex index;Row current;Throwable failure;boolean success;
	}
	private static final class Capture{Capture(int m,int a,int z,int s,int e,float id){model=m;ws=a;we=z;start=s;stop=e;identity=id;}final int model,ws,we,start,stop;final float identity;boolean overlap;}
	private static final class Row{String id,arm,host,split,source,status;byte[] bases;int start,stop,length;float gc;Capture chosen;final ArrayList<Capture> calls=new ArrayList<Capture>();}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final ArrayList<Row> rows=new ArrayList<Row>();private String ledger,fasta,prefix;private int expected;
}
