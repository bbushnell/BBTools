package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntList;

/** Joins final one-to-one false negatives to observed production stages.
 * This is evidence classification, not an exhaustive-model counterfactual:
 * partial shortlists remain explicitly unresolved until separately checked.
 * @author Raiden
 */
public final class Euk5sStageFailureSummary {
	public static void main(String[] args){
		String loci=null,stages=null,out=null,assembly=null,summaryIn=null;int expected=-1,models=-1,expectedAssemblies=-1;
		for(String arg : args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),value=arg.substring(eq+1);
			if(key.equals("loci")){loci=value;}else if(key.equals("stages")){stages=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("assembly")){assembly=value;}
			else if(key.equals("expectedfn")){expected=Integer.parseInt(value);}
			else if(key.equals("models")){models=Integer.parseInt(value);}
			else if(key.equals("summaryin")){summaryIn=value;}
			else if(key.equals("expectedassemblies")){expectedAssemblies=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+key);}
		}
		if(summaryIn!=null){
			require(out!=null && expected>=0 && expectedAssemblies>0,"Summary mode requires output and exact assembly/FN denominators");
			summarize(summaryIn,out,expectedAssemblies,expected);return;
		}
		require(loci!=null && stages!=null && out!=null && assembly!=null && expected>=0 && models>0,
			"Required loci, stages, outprefix, assembly, expectedfn and models");
		final Euk5sStageFailureSummary tool=new Euk5sStageFailureSummary(assembly,models);
		tool.loadLoci(loci,expected);tool.loadStages(stages);tool.write(out);
		System.err.println("EUK5S_STAGE_FAILURE_SUMMARY_PASS assembly="+assembly+" FN="+tool.count);
	}
	/** Pool existing classified rows; every assembly must contribute each category
	 * exactly once and its declared TOTAL must equal the category sum. */
	public static void summarize(String path,String prefix,int expectedAssemblies,long expectedFn){
		final ObjectMap<String,StageTotals> assemblies=new ObjectMap<String,StageTotals>(String.class,StageTotals.class);
		final ArrayList<StageTotals> rows=new ArrayList<StageTotals>();final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{Euk5sEndpointErrorSummary.header(in,"clade\tassembly\tstageEvidence\tFN",path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==4,"Summary input must retain clade, assembly, stage and FN count");
				final String clade=p.parseString(0),assembly=p.parseString(1);require(!clade.isEmpty() && !assembly.isEmpty(),"Stage totals require explicit cohort and assembly identities");
				StageTotals row=assemblies.get(assembly);if(row==null){row=new StageTotals(clade);assemblies.put(assembly,row);rows.add(row);}require(row.clade.equals(clade),"One assembly cannot be assigned to two clades");
				int category=-1;for(int i=0;i<CATEGORIES.length;i++){if(p.termEquals(CATEGORIES[i],2)){category=i;break;}}
				if(p.termEquals("TOTAL",2)){category=CATEGORIES.length;}require(category>=0 && !row.seen[category],"Unknown or duplicate per-assembly stage category");
				row.seen[category]=true;row.counts[category]=p.parseLong(3);require(row.counts[category]>=0,"Negative stage count cannot describe a partition");
			}
		}finally{require(!in.close(),"Stage summary read failed");}
		final ObjectMap<String,StageTotals> byClade=new ObjectMap<String,StageTotals>(String.class,StageTotals.class);final ArrayList<StageTotals> groups=new ArrayList<StageTotals>();long total=0;
		for(StageTotals row:rows){long sum=0;for(int i=0;i<row.counts.length;i++){require(row.seen[i],"Every assembly must provide zero categories explicitly as well as TOTAL");if(i<CATEGORIES.length){sum+=row.counts[i];}}
			require(sum==row.counts[CATEGORIES.length],"Stage categories disagree with their declared per-assembly FN total");total+=sum;
			StageTotals group=byClade.get(row.clade);if(group==null){group=new StageTotals(row.clade);byClade.put(row.clade,group);groups.add(group);}group.assemblies++;
			for(int i=0;i<row.counts.length;i++){group.counts[i]+=row.counts[i];}
		}
		require(rows.size()==expectedAssemblies && total==expectedFn,"Cohort summary must conserve the frozen assembly and FN denominators");
		RrnaCmWindowParity.fresh(prefix,".summary.tsv");final ByteStreamWriter writer=RrnaSourceTargetPrep.writer(prefix+".summary.tsv");final ByteBuilder b=new ByteBuilder();
		try{writer.println("clade\tassemblies\tstageEvidence\tFN");for(StageTotals group:groups){for(int i=0;i<group.counts.length;i++){
			writer.print(b.clear().append(group.clade).tab().append(group.assemblies).tab().append(i<CATEGORIES.length?CATEGORIES[i]:"TOTAL").tab().append(group.counts[i]).nl());
		}}}finally{RrnaCmWindowParity.finish(writer);}
		System.err.println("EUK5S_STAGE_COHORT_SUMMARY_PASS assemblies="+rows.size()+" FN="+total);
	}
	private static final class StageTotals{
		StageTotals(String clade_){clade=clade_;}
		final String clade;int assemblies;final long[] counts=new long[CATEGORIES.length+1];final boolean[] seen=new boolean[CATEGORIES.length+1];
	}
	public Euk5sStageFailureSummary(String assembly_,int models_){
		this(assembly_,models_,"euk5S");
	}
	Euk5sStageFailureSummary(String assembly_,int models_,String family_){
		require(assembly_!=null && !assembly_.isEmpty() && models_>0,"An exact assembly and model count bind the diagnostic");
		require(family_!=null && !family_.isEmpty(),"Stage evidence must bind the caller family explicitly");
		assembly=assembly_;models=models_;family=family_;
	}
	public void loadLoci(String path,int expected){
		require(groups.isEmpty() && count==0,"Load one assembly once");
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> unique=new ObjectSet<String>(String.class);
		try{
			Euk5sEndpointErrorSummary.header(in,LOCUS_HEADER,path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==12 && p.termEquals(assembly,0),"Locus schema/assembly mismatch");
				if(!p.termEquals("TRUTH",2) || !p.termEquals("FN",7)){continue;}
				require(p.parseInt(8)==-1 && p.termEquals("NA",9) && p.termEquals("NA",10),"False negative cannot have a matched partner");
				final String seq=p.parseString(1);final int start=p.parseInt(4)-1,stop=p.parseInt(5)-1;
				final int strand=p.termEquals("+",6) ? 0 : p.termEquals("-",6) ? 1 : -1;
				require(strand>=0 && start>=0 && stop>=start,"False negative needs valid genomic coordinates");
				require(unique.add(seq+"\t"+strand+"\t"+start+"\t"+stop),"Duplicate FN coordinates");
				Group group=bySeq.get(seq);
				if(group==null){group=new Group(seq);bySeq.put(seq,group);groups.add(group);}
				group.pending.add(new Evidence(start,stop,strand));count++;
			}
		}finally{require(!in.close(),"FN loci read failed");}
		require(count==expected,"FN row count differs from the accepted grade: "+count+" != "+expected);
		for(Group group : groups){group.finish();}
	}
	public void loadStages(String path){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		final IntList seeds=new IntList();String lastName=null;Group group=null;
		try{
			Euk5sEndpointErrorSummary.header(in,NcrnaStageDiagLogSink.HEADER,path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==21,"Stage diagnostic must retain21 columns");
				if(!p.termEquals(family,1)){continue;}
				// Per-model identities are separate evidence, not terminal windows or
				// extra candidates. Keep existing FN stage counts unchanged.
				if(p.termEquals("MODEL",0)){continue;}
				if(lastName==null || !p.termEquals(lastName,2)){
					lastName=p.parseString(2);group=bySeq.get(RrnaCallerModelAssignment.first(lastName));
				}
				if(group==null){continue;}
				final int strand=p.parseInt(3),length=p.parseInt(4);
				require((strand==0 || strand==1) && length>0,"Stage strand/contig length invalid");
				if(group.length<0){group.length=length;}else{require(group.length==length,"Conflicting contig lengths");}
				if(p.termEquals("SEED",0)){
					require(!group.seedSeen[strand],"Duplicate seed stream for one contig/strand");group.seedSeen[strand]=true;
					parseSeeds(row,p,seeds,length);
					for(Evidence e : group.loci){if(e.strand==strand){
						require(e.stop<length,"FN extends beyond diagnostic contig");
						final int first=(strand==0 ? e.start : length-1-e.stop)+8;
						final int last=(strand==0 ? e.stop : length-1-e.start)-8;
						e.localSeeds=last<first ? 0 : lowerBound(seeds,last+1)-lowerBound(seeds,first);
					}}
				}else if(p.termEquals("WINDOW",0)){
					final int a=p.parseInt(6),b=p.parseInt(7);validSpan(a,b,length);
					final int start=strand==0 ? a : length-1-b,stop=strand==0 ? b : length-1-a;
					final String outcome=p.parseString(8);final int shortlist=p.parseInt(10);
					require(shortlist>=-1 && shortlist<=models,"Shortlist exceeds the bound library");
					final int ca=p.parseInt(16),cb=p.parseInt(17);
					final boolean accepted=outcome.startsWith("ACCEPT_");
					if(accepted){validSpan(ca,cb,length);}
					final int callStart=strand==0 ? ca : length-1-cb,callStop=strand==0 ? cb : length-1-ca;
					for(int i=group.firstOverlap(start);i<group.loci.length && group.loci[i].start<=stop;i++){
						final Evidence e=group.loci[i];if(e.strand!=strand || e.stop<start){continue;}
						e.overlapping++;
						if(start<=e.start && stop>=e.stop){
							e.containing++;
							if(accepted){e.containingAccepts++;}
							else if(outcome.equals("REJECT_KHITS")){e.seedRejected++;}
							else if(outcome.equals("REJECT_INDEX")){e.indexRejected++;}
							else if(outcome.equals("REJECT_VERIFY") || outcome.equals("REJECT_ID") || outcome.equals("REJECT_RESCUE")){
								if(shortlist==models){e.allModelRejected++;}else{e.partialModelRejected++;}
							}else if(outcome.equals("REJECT_LEN") || outcome.equals("REJECT_MAXLEN") || outcome.equals("REJECT_WLEN")){e.lengthRejected++;}
							else{throw new IllegalArgumentException("Unknown terminal window outcome: "+outcome);}
						}
						if(accepted && reciprocal(e,callStart,callStop)){e.acceptedMatches++;}
					}
				}else if(p.termEquals("PATH",0)){
					final int a=p.parseInt(16),b=p.parseInt(17);validSpan(a,b,length);
					final int start=strand==0 ? a : length-1-b,stop=strand==0 ? b : length-1-a;
					final boolean kept=p.termEquals("DP_RETAINED",8);
					require(kept || p.termEquals("DP_INPUT",8),"Unknown DP observation");
					for(int i=group.firstOverlap(start);i<group.loci.length && group.loci[i].start<=stop;i++){
						final Evidence e=group.loci[i];if(e.strand!=strand || !reciprocal(e,start,stop)){continue;}
						if(kept){e.retainedMatches++;}else{e.dpMatches++;}
					}
				}else{throw new IllegalArgumentException("Unknown stage event");}
			}
		}finally{require(!in.close(),"Stage observations read failed");}
	}
	private static void parseSeeds(byte[] row,LineParser1 p,IntList seeds,int length){
		seeds.clear();final int expected=p.parseInt(19);require(expected>=0,"Negative seed count");
		if(expected==0){require(p.termEquals(".",20),"Zero seeds need explicit empty field");return;}
		p.setBounds(20);final int end=p.b();int value=0,digits=0;
		for(int i=p.a();i<=end;i++){
			if(i==end || row[i]==','){
				require(digits>0 && value>=8 && value<=length-9,"Invalid complete17mer center");
				require(seeds.size==0 || seeds.array[seeds.size-1]<value,"Seed centers must be strictly increasing");
				seeds.add(value);value=0;digits=0;
			}else{final int d=row[i]-'0';require(d>=0 && d<=9 && value<=(Integer.MAX_VALUE-d)/10,"Malformed seed position");value=value*10+d;digits++;}
		}
		require(seeds.size==expected,"Seed count differs from its positions");
	}
	private static int lowerBound(IntList values,int key){
		assert(values!=null) : "Binary search uses the validated sorted seed stream";
		int lo=0,hi=values.size;while(lo<hi){final int mid=(lo+hi)>>>1;if(values.array[mid]<key){lo=mid+1;}else{hi=mid;}}return lo;
	}
	static boolean reciprocal(Evidence e,int start,int stop){
		assert(e!=null && stop>=start) : "Candidate overlap requires an actual span";
		final long overlap=Math.max(0,Math.min(e.stop,stop)-Math.max(e.start,start)+1);
		return 2*overlap>=e.stop-e.start+1L && 2*overlap>=stop-start+1L;
	}
	public static String category(Evidence e){
		assert(e!=null) : "Classify one final unmatched truth locus";
		if(e.retainedMatches>0){return "MATCHING_COMPETITION";}
		if(e.dpMatches>0){return "DP_DROPPED_QUALIFYING_CANDIDATE";}
		if(e.acceptedMatches>0){return "PRE_DP_RESOLUTION";}
		if(e.containingAccepts>0){return "ACCEPTED_OTHER_LOCUS";}
		if(e.partialModelRejected>0){return "SHORTLIST_OR_IDENTITY_UNRESOLVED";}
		if(e.allModelRejected>0){return "ALL_MODELS_VERIFICATION_REJECTED";}
		if(e.lengthRejected>0){return "LENGTH_REJECTED";}
		if(e.indexRejected>0){return "SHORTLIST_REJECTED";}
		if(e.seedRejected>0){return "WINDOW_SEED_GATE_REJECTED";}
		if(e.containing>0){throw new IllegalStateException("Containing window lacks a classified terminal outcome");}
		if(e.overlapping>0){return "WINDOW_TRUNCATED";}
		if(e.localSeeds<0){return "NO_SEED_OBSERVATION";}
		return e.localSeeds<2 ? "FEWER_THAN_TWO_LOCAL_SEEDS_NO_WINDOW" : "WINDOW_ABSENT";
	}
	public void write(String prefix){
		RrnaCmWindowParity.fresh(prefix,new String[]{".loci.tsv.gz",".summary.tsv"});
		final ByteStreamWriter loci=RrnaSourceTargetPrep.writer(prefix+".loci.tsv.gz"),summary=RrnaSourceTargetPrep.writer(prefix+".summary.tsv");
		final long[] counts=new long[CATEGORIES.length];final ByteBuilder b=new ByteBuilder();long written=0;
		try{
			loci.println("assembly\tseqid\tstart\tstop\tstrand\tlocal17merHits\toverlappingWindows\tcontainingWindows\tacceptedMatches\tdpInputMatches\tdpRetainedMatches\tcontainingAccepts\tallModelRejects\tpartialModelRejects\tindexRejects\tseedRejects\tlengthRejects\tstageEvidence");
			for(Group group : groups){for(Evidence e : group.loci){
				require(e.retainedMatches<=e.dpMatches,"Retained matching events exceed actual DP input events");
				final String category=category(e);int index=-1;for(int i=0;i<CATEGORIES.length;i++){if(CATEGORIES[i].equals(category)){index=i;break;}}
				require(index>=0,"Unregistered category");counts[index]++;written++;
				b.clear().append(assembly).tab().append(group.name).tab().append(e.start+1).tab().append(e.stop+1).tab().append(e.strand==0?'+':'-').tab()
					.append(e.localSeeds).tab().append(e.overlapping).tab().append(e.containing).tab().append(e.acceptedMatches).tab().append(e.dpMatches).tab()
					.append(e.retainedMatches).tab().append(e.containingAccepts).tab().append(e.allModelRejected).tab().append(e.partialModelRejected).tab()
					.append(e.indexRejected).tab().append(e.seedRejected).tab().append(e.lengthRejected).tab().append(category).nl();loci.print(b);
			}}
			require(written==count,"Stage categories must conserve the final FN denominator");
			summary.println("assembly\tstageEvidence\tFN");
			for(int i=0;i<CATEGORIES.length;i++){summary.print(b.clear().append(assembly).tab().append(CATEGORIES[i]).tab().append(counts[i]).nl());}
			summary.print(b.clear().append(assembly).append("\tTOTAL\t").append(count).nl());
		}finally{RrnaCmWindowParity.finish(loci,summary);}
	}
	static void validSpan(int a,int b,int length){require(a>=0 && b>=a && b<length,"Diagnostic span lies outside the contig");}
	static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	public static final class Evidence{
		public Evidence(int a,int b,int s){start=a;stop=b;strand=s;}
		final int start,stop,strand;public int localSeeds=-1;
		public int overlapping,containing,acceptedMatches,dpMatches,retainedMatches,containingAccepts;
		int allModelRejected,partialModelRejected,indexRejected,seedRejected,lengthRejected;
	}
	private static final class Group{
		Group(String name_){name=name_;}
		void finish(){
			assert(!pending.isEmpty()) : "Only contigs with final false negatives enter the join";
			loci=pending.toArray(new Evidence[pending.size()]);pending.clear();
			Arrays.sort(loci,(a,b)->Integer.compare(a.start,b.start));for(Evidence e : loci){maxLength=Math.max(maxLength,e.stop-e.start+1);}
		}
		int firstOverlap(int start){
			assert(loci!=null && maxLength>0) : "Build the sorted interval index before reading traces";
			final long key=(long)start-maxLength+1;int lo=0,hi=loci.length;
			while(lo<hi){final int mid=(lo+hi)>>>1;if(loci[mid].start<key){lo=mid+1;}else{hi=mid;}}return lo;
		}
		final String name;final ArrayList<Evidence> pending=new ArrayList<Evidence>();
		final boolean[] seedSeen=new boolean[2];Evidence[] loci;int maxLength,length=-1;
	}
	final String assembly,family;final int models;long count;
	private final ObjectMap<String,Group> bySeq=new ObjectMap<String,Group>(String.class,Group.class);
	private final ArrayList<Group> groups=new ArrayList<Group>();
	private static final String LOCUS_HEADER="assembly\tseqid\trole\tindex\tstart\tstop\tstrand\toutcome\tpartnerIndex\td5\td3\tcontrolOverlapAnyStrand";
	private static final String[] CATEGORIES={"MATCHING_COMPETITION","DP_DROPPED_QUALIFYING_CANDIDATE","PRE_DP_RESOLUTION",
		"ACCEPTED_OTHER_LOCUS","SHORTLIST_OR_IDENTITY_UNRESOLVED","ALL_MODELS_VERIFICATION_REJECTED","LENGTH_REJECTED",
		"SHORTLIST_REJECTED","WINDOW_SEED_GATE_REJECTED","WINDOW_TRUNCATED","NO_SEED_OBSERVATION",
		"FEWER_THAN_TWO_LOCAL_SEEDS_NO_WINDOW","WINDOW_ABSENT"};
}
