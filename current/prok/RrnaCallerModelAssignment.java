package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.LongHashSet;
import map.ObjectSet;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import structures.ByteBuilder;
import structures.ListNum;

/** Offline pre-DP model dispatch on independently labeled, sense-oriented genes.
 * Reuses production family controls and scavenger; never substitutes taxonomy or
 * exhaustive best identity for a failed seed/window/shortlist/verifier path.
 * Launchers bind explicit consensus and seed input files by canonical sha80.
 * @author Ganyu
 */
public final class RrnaCallerModelAssignment {

	public static void main(String[] args){
		String in=null,targets=null,exact=null,consensus=null,seeds=null,prefix=null;int expected=-1,models=-1,flank=-1;boolean requireProductionSeeds=false;
		for(String arg:args){final int eq=arg.indexOf('=');require(eq>0,"Expected key=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("targets")){targets=v;}else if(k.equals("exact")){exact=v;}else if(k.equals("consensus")){consensus=v;}else if(k.equals("seeds")){seeds=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("expectedmodels")){models=Integer.parseInt(v);}else if(k.equals("flank")){flank=Integer.parseInt(v);}
			else if(k.equals("requireproductionseeds")){require(v.equals("t")||v.equals("f"),"requireproductionseeds must be t or f");requireProductionSeeds=v.equals("t");}else{throw new IllegalArgumentException(k);}}
		require(in!=null&&targets!=null&&exact!=null&&consensus!=null&&seeds!=null&&prefix!=null&&expected>0&&models>0&&flank>0,"Explicit input/target/exact/library/seeds/output/count/flank contract required");
		require(GeneCaller.ncrnaFamilies.isEmpty(),"Dedicated offline process must begin with no loaded families");
		CallGenes.loadRrna17Resources("euk5s");require(GeneCaller.ncrnaFamilies.size()==1,"Exactly one production euk5S family required");
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensus);final String[] names=TrnaConsensusBuilder.lastLoadedNames.clone();
		require(library.length==models&&names.length==models,"Explicit consensus count changed");
		final LongHashSet seedSet=ProkObject.loadLongKmers(seeds,17);
		if(requireProductionSeeds){require(sameSeeds(seedSet,GeneCaller.ncrnaFamilies.get(0).kmerSet),"Explicit arm-A seeds differ from effective production seed membership");}
		new RrnaCallerModelAssignment(GeneCaller.ncrnaFamilies.get(0),library,names,seedSet).run(in,targets,exact,prefix,expected,flank);
	}

	public RrnaCallerModelAssignment(NcrnaFamily f,byte[][] library,String[] labels,LongHashSet seeds){
		require(f.name.equals("euk5S")&&f.kLong==17&&f.models==null&&f.boundary5NetTemplate==null&&f.boundary3NetTemplate==null&&f.boundaryNetsByModel==null,
			"This offline adapter supports the current NN-off/HBM-free euk5S family only");
		require(library.length>0&&labels.length==library.length,"Every candidate consensus needs a name");
		require(seeds!=null&&seeds.size()>0,"The explicit candidate seed set must be nonempty");seedCount=seeds.size();
		productionSeedCount=f.kmerSet.size();productionSeedsEqual=sameSeeds(seeds,f.kmerSet);
		System.err.println("ASSIGNMENT_SEEDS explicit="+seedCount+" production="+productionSeedCount+" sameMembership="+productionSeedsEqual);
		names=new String[labels.length];final ObjectSet<String> unique=new ObjectSet<String>(String.class);
		for(int i=0;i<labels.length;i++){names[i]=first(labels[i]);require(unique.add(names[i])&&library[i]!=null&&library[i].length>0,"Unique nonempty models required");}
		family=f;
		scavenger=new NcrnaScavenger(library,null,names,seeds,f.kLong,f.minLen,f.windowPad,
			f.indexK,f.indexTopN,f.adaptive,f.adaptFloor,f.adaptTopFrac,f.adaptQFrac,f.fixedMinHits,
			f.scoreA,f.scoreB,f.idPass,f.idBorderline,f.boundary5NetTemplate,f.boundary3NetTemplate,
			f.boundaryStartTable,f.boundaryStopTable,f.boundaryStartInside,f.boundaryStartOutside,f.boundaryStopInside,f.boundaryStopOutside,
			f.boundaryMeanLen,f.boundaryStartOffsets,f.boundaryStopOffsets,f.boundaryMarginStart,f.boundaryMarginStop);
		GeneCaller.applyNcrnaFamilyControls(scavenger,f);
		index=new ConservedRnaSeedIndex(new LongHashSet[]{seeds},new boolean[]{f.voteWindows||f.voteEnds});
		scavenger.setInstrumentSink((contig,strand,model,ws,we,start,stop,copy,offset,trim,nn,source,identity,length)->{
			require(!nn,"No endpoint network may run while generating training assignments");
			captures.add(new Capture(strand,model,ws,we,start,stop,identity));
		});
	}

	public Result assign(Read r,int left,int right){
		assert(r!=null&&r.bases!=null):"Assignment operates on a real flanked gene";
		require(left>=0&&right>=0&&left+right<r.length(),"Invalid actual flanks");
		captures.clear();final int[] hitCounts=new int[2];final long before=scavenger.alignmentCount();int calls=0;
		for(int strand=0;strand<2;strand++){
			final byte[] bases=strand==0?r.bases:AminoAcid.reverseComplementBases(r.bases);
			final ConservedRnaSeedIndex.ScanResult hits=index.scan(bases);hitCounts[strand]=hits.hits(0).length;
			// Same bit-index hit stream and family-local call used by GeneCaller.makeRnas.
			final ArrayList<Orf> found=scavenger.scavenge(first(r.id),bases,strand,new ArrayList<int[]>(),hits.hits(0),hits.keys(0));
			calls+=found.size();
		}
		require(calls==captures.size(),"Instrumentation must conserve the returned pre-DP calls");
		return select(new ArrayList<Capture>(captures),left,r.length()-right-1,hitCounts,scavenger.alignmentCount()-before);
	}

	public static Result select(ArrayList<Capture> all,int start,int stop,int[] hits,long alignments){
		assert(start>=0&&stop>=start&&hits.length==2):"Truth is in the sense-oriented record frame";
		Capture chosen=null;int eligible=0;
		for(Capture c:all){final int overlap=Math.max(0,Math.min(stop,c.stop)-Math.max(start,c.start)+1);
			if(c.strand==0&&2L*overlap>=Math.min(stop-start+1,c.stop-c.start+1)){eligible++;chosen=c;}}
		final String status=eligible>1?"AMBIGUOUS":eligible==1?"ASSIGNED":all.isEmpty()?"NO_ACCEPTED_CALL":"NO_MATCHING_CALL";
		return new Result(status,eligible==1?chosen:null,all,hits,alignments);
	}

	public void run(String in,String targetPath,String exact,String prefix,int expected,int flank){
		for(String suffix:new String[]{".assignments.tsv",".candidates.tsv.gz",".assigned.fa.gz",".summary.tsv",".config.tsv"}){require(!new File(prefix+suffix).exists(),"Output already exists");}
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		final RrnaTargetExtractionCheck truth=new RrnaTargetExtractionCheck();truth.loadTargets(targetPath);truth.loadExact(exact,"5S");
		require(truth.byId.size()==expected,"Target count changed");
		final ByteStreamWriter rows=writer(prefix+".assignments.tsv"),events=writer(prefix+".candidates.tsv.gz"),fasta=writer(prefix+".assigned.fa.gz");
		rows.println("target_id\tassembly\tstatus\tmodel\tidentity\twindowStart0\twindowStop0\tcallStart0\tcallStop0\tleftEdgeContact\trightEdgeContact\tleftArtificialEdge\trightArtificialEdge\tforwardSeedHits\treverseSeedHits\talignments\tacceptedCalls");
		events.println("target_id\tstrand\tmodel\tidentity\twindowStart0\twindowStop0\tcallStart0\tcallStop0\tleftEdgeContact\trightEdgeContact");
		final Streamer input=RrnaTargetExtractionCheck.stream(in);input.start();final ByteBuilder b=new ByteBuilder();int total=0,assigned=0,ambiguous=0,noCall=0,unmatched=0,edge=0,artificial=0;boolean finished=false;
		try{for(ListNum<Read> batch;(batch=input.nextList())!=null;){for(Read r:batch.list){
			final String id=first(r.id);final int pipe=id.indexOf('|',4);require(id.startsWith("tid|")&&pipe>4,"Checked FASTA requires attributed target IDs");
			final String targetId=id.substring(pipe+1);final RrnaTargetExtractionCheck.Target t=truth.byId.get(targetId);
			require(t!=null&&!t.seen&&t.family.equals("5S")&&t.exact!=null&&Integer.parseInt(id.substring(4,pipe))==t.taxid,"Foreign/duplicate/non5S/missing exact target");t.seen=true;
			require(field(r.id,"assembly_id").equals(t.assembly)&&field(r.id,"source_contig").equals(t.seq)&&field(r.id,"source_gff").equals(t.source)
				&&number(r.id,"source_start")==t.start&&number(r.id,"source_end")==t.stop&&field(r.id,"source_strand").equals(t.strand),"Source attribution changed");
			final int left=number(r.id,"lflank"),right=number(r.id,"rflank");require(left<=flank&&right<=flank,"Actual flank exceeds frozen request");
			NcrnaAnchorSlop.requireExactGene(r.bases,left,right,t.exact);
			final Result result=assign(r,left,right);total++;
			if(result.status.equals("ASSIGNED")){assigned++;}else if(result.status.equals("AMBIGUOUS")){ambiguous++;}else if(result.status.equals("NO_ACCEPTED_CALL")){noCall++;}else{unmatched++;}
			b.clear().append(targetId).tab().append(t.assembly).tab().append(result.status);
			final Capture c=result.chosen;
			if(c==null){b.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}
			else{final boolean le=c.ws==0,re=c.we==r.length()-1,la=le&&left==flank,ra=re&&right==flank;
				if(le||re){edge++;}if(la||ra){artificial++;}
				b.tab().append(names[c.model]).tab().append(c.identity,8).tab().append(c.ws).tab().append(c.we).tab().append(c.start).tab().append(c.stop).tab().append(le?1:0).tab().append(re?1:0).tab().append(la?1:0).tab().append(ra?1:0);
				fasta.print(new ByteBuilder().append('>').append(r.id).append(" accepted_model=").append(names[c.model]).append(" model_assignment_source=production_scavenger_pre_dp").nl().append(r.bases).nl());}
			b.tab().append(result.hits[0]).tab().append(result.hits[1]).tab().append(result.alignments).tab().append(result.all.size()).nl();rows.print(b);
			for(Capture e:result.all){b.clear().append(targetId).tab().append(e.strand).tab().append(names[e.model]).tab().append(e.identity,8).tab().append(e.ws).tab().append(e.we).tab().append(e.start).tab().append(e.stop).tab().append(e.ws==0?1:0).tab().append(e.we==r.length()-1?1:0).nl();events.print(b);}
		}input.returnList(batch);}require(!input.errorState()&&total==expected,"Input failure or target conservation failure");finished=true;
		}finally{if(!finished){input.close();}RrnaCmWindowParity.finish(rows,events,fasta);}
		final ByteStreamWriter summary=writer(prefix+".summary.tsv");
		try{summary.println("targets\tassigned\tambiguous\tnoAcceptedCall\tnoMatchingCall\tassignedEdgeContact\tassignedArtificialEdgeContact");summary.print(b.clear().append(total).tab().append(assigned).tab().append(ambiguous).tab().append(noCall).tab().append(unmatched).tab().append(edge).tab().append(artificial).nl());}finally{RrnaCmWindowParity.finish(summary);}
		final ByteStreamWriter config=writer(prefix+".config.tsv");
		try{config.println("family\tseedK\tseedMinHits\tmapK\ttopN\tmapMinHits\tidPass\twindowPad\tminLen\tcollapseFrac\tpass2\treuseQuantum\tmodels\texplicitSeedCount\tproductionSeedCount\tsameProductionMembership");
			config.println(family.name+"\t"+family.kLong+"\t"+family.seedMinHits+"\t"+family.indexK+"\t"+family.indexTopN+"\t"+family.fixedMinHits+"\t"+family.idPass+"\t"+family.windowPad+"\t"+family.minLen+"\t"+family.collapseFrac+"\t"+family.scavengePass2+"\t"+family.reuseConsensusAlignment+"\t"+names.length+"\t"+seedCount+"\t"+productionSeedCount+"\t"+productionSeedsEqual);
		}finally{RrnaCmWindowParity.finish(config);}
		System.err.println("CALLER_MODEL_ASSIGNMENT_PASS targets="+total+" assigned="+assigned+" artificialEdgeContact="+artificial+"; pre-DP, flank-limited context; no genomic-recall claim");
	}

	static String field(String h,String k){return RrnaTargetExtractionCheck.field(h,k);}
	public static boolean sameSeeds(LongHashSet a,LongHashSet b){require(a!=null&&b!=null,"Seed comparison requires both actual sets");if(a.size()!=b.size()){return false;}for(long key:a.toArray()){if(!b.contains(key)){return false;}}return true;}
	static int number(String h,String k){return Integer.parseInt(field(h,k));}
	static String first(String s){require(s!=null&&!s.isEmpty(),"Nonempty name required");int n=0;while(n<s.length()&&!Character.isWhitespace(s.charAt(n))){n++;}return s.substring(0,n);}
	static ByteStreamWriter writer(String path){return RrnaCmWindowParity.writer(path);}
	static void require(boolean b,String why){if(!b){throw new IllegalArgumentException(why);}}
	public static final class Capture{public Capture(int strand_,int model_,int ws_,int we_,int start_,int stop_,float identity_){require(ws_>=0&&we_>=ws_&&start_>=ws_&&stop_<=we_&&stop_>=start_,"Invalid production-capture bounds");strand=strand_;model=model_;ws=ws_;we=we_;start=start_;stop=stop_;identity=identity_;}public final int strand,model,ws,we,start,stop;final float identity;}
	public static final class Result{Result(String s,Capture c,ArrayList<Capture> a,int[] h,long n){status=s;chosen=c;all=a;hits=h;alignments=n;}public final String status;public final Capture chosen;final ArrayList<Capture> all;final int[] hits;final long alignments;}
	final int seedCount,productionSeedCount;final boolean productionSeedsEqual;final NcrnaFamily family;final String[] names;final NcrnaScavenger scavenger;final ConservedRnaSeedIndex index;final ArrayList<Capture> captures=new ArrayList<Capture>();
}
