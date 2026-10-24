package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Reports saved tier counters without calling genes or changing grading.
 * Missing counters are unavailable, not zero. An explicit original-only baseline
 * binding alone permits deriving its original work from the family16S total.
 * @author Brian, Ganyu */
public final class NcrnaMitoWorkReport {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoWorkReport.class,false);
		try{final NcrnaMitoWorkReport x=new NcrnaMitoWorkReport();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] p=Parse.splitOnFirst(arg,'=');require(p.length==2,"Use flag=value");final String key=p[0].toLowerCase(Locale.ROOT);
			require(KNOWN.contains("|"+key+"|") && flags.put(key,p[1])==null,"Unknown/repeated argument: "+key);}
		expected=Parse.parseIntKMG(flags.get("expectedruns")==null?"8":flags.get("expectedruns"));require(expected>0 && (expected&1)==0,"Require a positive even count of paired baseline/C runs");
		pins=EukPanelGffValidation.hashes(required("verified"));bound(required("runs"));prefix=required("outprefix");RrnaResourceIO.fresh(prefix,".tiers.tsv",".inputs.tsv");Shared.setThreads(1);ReadWrite.setZipThreads(1);
	}
	private void run(){
		final ByteFile in=RrnaResourceIO.open(required("runs"));final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class),roots=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,ROSTER_HEADER,required("runs"));for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==5,"Tier roster requires arm/panel/result_root/seal_sha80/counter_policy");
			final Run r=new Run(p);require(seen.add(r.panel+'\t'+r.arm) && roots.add(r.root),"Run identities and source roots must be unique");
			require(bound(r.root+"/MANIFEST.sha80").equals(r.seal),"Run manifest must match the explicitly bound roster seal");
			final ObjectMap<String,String> manifest=EukPanelGffValidation.hashes(r.root+"/MANIFEST.sha80");r.logHash=bound(r.root+"/stderr");require(r.logHash.equals(manifest.get("stderr")),"Caller stderr must be the payload bound by its own run manifest");
			readCounters(r);runs.add(r);
		}}finally{require(!in.close(),"Tier roster read failed");}
		require(runs.size()==expected,"Saved caller run count differs from expectedruns");for(Run r:runs){require(seen.contains(r.panel+'\t'+(r.arm.equals("baseline")?"C":"baseline")),"Every panel needs both explicitly bound arms");}
		write();
	}
	private static void readCounters(Run r){
		final ByteFile in=RrnaResourceIO.open(r.root+"/stderr");
		try{for(byte[] line;(line=in.nextLine())!=null;){for(int i=0;i<PREFIXES.length;i++){if(Tools.startsWith(line,PREFIXES[i])){
				require(!r.seen[i],"Repeated tier/family counter: "+PREFIXES[i]);r.seen[i]=true;r.counts[i]=counter(line,PREFIXES[i].length());break;
			}}}}finally{require(!in.close(),"Caller stderr read failed: "+r.root);}
		final int tierFields=(r.seen[1]?1:0)+(r.seen[2]?1:0)+(r.seen[3]?1:0);require(tierFields==0 || tierFields==3,"A partially emitted tier-counter trio is corrupt, not a zero count");
		if(tierFields==3){
			require(r.seen[0] && r.counts[1]<=Long.MAX_VALUE-r.counts[2] && r.counts[1]+r.counts[2]==r.counts[0],"Original plus mitochondrial alignments must equal the measured family16S total");
			require(!r.baseline || r.counts[2]==0 && r.counts[3]==0,"An ORIGINAL_ONLY baseline cannot report mitochondrial work or skips");r.state="MEASURED";
		}else if(r.baseline && r.seen[0]){r.counts[1]=r.counts[0];r.counts[2]=r.counts[3]=0;r.state="BASELINE_DERIVED";}
		else{r.state="UNAVAILABLE";}
	}
	/** Exact unsigned decimal parse; overflow and nonnumeric diagnostics fail loud. */
	private static long counter(byte[] line,int from){
		int end=line.length;while(from<end && line[from]<=32){from++;}while(end>from && line[end-1]<=32){end--;}require(from<end,"Missing counter value");long n=0;
		for(int i=from;i<end;i++){final int digit=line[i]-'0';require(digit>=0 && digit<=9 && n<=(Long.MAX_VALUE-digit)/10,"Counter must be a nonnegative long integer");n=n*10+digit;}return n;
	}
	private void write(){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".tiers.tsv"),inputs=RrnaResourceIO.writer(prefix+".inputs.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("arm\tpanel\tcounter_policy\tfamily16S_alignments\toriginal16S_alignments\tmito_alignments\tmito_prealignment_length_skips\tstate\tinterpretation");inputs.println("role\tpath\tsha80");
			inputs.print(b.clear().append("runs\t").append(required("runs")).tab().append(bound(required("runs"))).nl());
			for(Run r:runs){b.clear().append(r.arm).tab().append(r.panel).tab().append(r.policy);for(int i=0;i<4;i++){b.tab();if(i==0?r.seen[0]:!r.state.equals("UNAVAILABLE")){b.append(r.counts[i]);}else{b.append("NA");}}
				out.print(b.tab().append(r.state).append("\tFINAL_CALLER_WORK_NOT_BIOLOGICAL_COUNTS\n"));
				inputs.print(b.clear().append("manifest\t").append(r.root+"/MANIFEST.sha80").tab().append(r.seal).nl());inputs.print(b.clear().append("stderr\t").append(r.root+"/stderr").tab().append(r.logHash).nl());}
		}finally{RrnaResourceIO.finish(out,inputs);}
		System.err.println("MITO_WORK_REPORT_PASS runs="+runs.size()+" new_caller_runs=0");
	}
	private String required(String key){final String v=flags.get(key);require(v!=null && !v.isEmpty(),"Missing argument: "+key);return v;}
	private String bound(String path){final String hash=pins.get(path);require(hash!=null,"Missing producer-verified input: "+path);return hash;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Run{
		Run(LineParser1 p){arm=p.parseString(0);panel=p.parseString(1);root=p.parseString(2);seal=p.parseString(3);policy=p.parseString(4);baseline=arm.equals("baseline");
			require((baseline || arm.equals("C")) && panel.matches("P[1-4]") && policy.equals(baseline?"ORIGINAL_ONLY":"COVERAGE"),"Baseline derivation requires explicit baseline/ORIGINAL_ONLY binding; candidate requires C/COVERAGE");
			require(root.startsWith("/") && !root.contains("/../") && !root.endsWith("/..") && SeedOffsetTable.validSha80(seal),"Each saved run requires an absolute bound root and sha80 seal");}
		final String arm,panel,root,seal,policy;final boolean baseline;final long[] counts=new long[4];final boolean[] seen=new boolean[4];String state,logHash;
	}
	public static final String ROSTER_HEADER="arm\tpanel\tresult_root\tseal_sha80\tcounter_policy";
	private static final String[] PREFIXES={"ncRNA Alignments (16S):","ncRNA Original16S Alignments:","ncRNA Mito Alignments:","ncRNA Mito Prealignment Length Skips:"};
	private static final String KNOWN="|runs|verified|outprefix|expectedruns|";
	private final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);private ObjectMap<String,String> pins;private String prefix;private int expected;private final ArrayList<Run> runs=new ArrayList<Run>();
}
