package prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Raw equal-pass/borderline caller-stage accounting over every declared input record.
 * MODEL observations include primary and physical-edge retry phases. WINDOW
 * shortlist size describes its terminal phase, not total entries across retries.
 * Events and distinct input records are separate. HBM evaluations are unobserved;
 * equal thresholds disable rescue, but edge verification may still score HBMs.
 * absent SEED observations are not silently converted into observed zero-hit events.
 * @author Raiden
 */
public final class NcrnaShredStages {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaShredStages.class, false);final NcrnaShredStages x=new NcrnaShredStages();
		String fasta=null, stages=null, gff=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){fasta=v;}else if(k.equals("stages")){stages=v;}else if(k.equals("gff")){gff=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("family")){x.family=v;}else if(k.equals("prefix")){x.prefix=v;}
			else if(k.equals("label")){x.label=v;}else if(k.equals("idpass")){x.cut=Float.parseFloat(v);}else if(k.equals("idborderline")){x.border=Float.parseFloat(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(fasta!=null && stages!=null && gff!=null && out!=null && expected>0 && x.family!=null && x.prefix!=null && x.label!=null
			&& Float.isFinite(x.cut) && x.cut>0 && x.cut<=1 && x.border==x.cut, "Bind all inputs and explicit equal raw pass/borderline thresholds; this reducer does not infer HBM evaluation");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{x.load(fasta, expected);x.stages(stages);x.calls(gff);x.write(out);}finally{Shared.closeStream(pp.outstream);}
	}
	private void load(String path, int expected){
		RrnaResourceIO.stream(path, r->{
			final String id=RrnaResourceIO.first(r.id);require(r.bases!=null && r.length()>0, "Every input record needs its full positive base denominator");
			final Member m=new Member(id, r.length());require(byId.put(id, m)==null, "Duplicate input identity: "+id);members.add(m);bases+=r.length();
		});require(members.size()==expected, "All declared input records must survive stage accounting");
	}
	private void stages(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, NcrnaStageDiagLogSink.HEADER, path);for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);require(p.terms()==21, "Caller stage schema changed");if(!p.termEquals(family, 1)){continue;}
			final Member m=byId.get(RrnaResourceIO.first(p.parseString(2)));require(m!=null && p.parseInt(4)==m.length, "Stage event must bind the exact input record and length");
			final int strand=p.parseInt(3);require(strand==0 || strand==1, "Stage strand must be explicit");
			if(p.termEquals("SEED", 0)){
				final int n=p.parseInt(19);require(n>=0, "Observed seed occurrence count cannot be negative");m.events[0]++;m.events[1]+=n;continue;
			}
			if(p.termEquals("PATH", 0)){
				if(p.termEquals("DP_INPUT", 8)){m.events[7]++;}else{require(p.termEquals("DP_RETAINED", 8), "Unknown path outcome");m.events[8]++;}continue;
			}
			final int pass=p.parseInt(5);require(pass==1 || pass==2, "Alignment events must identify actual caller pass1 or2");
			if(p.termEquals("MODEL", 0)){
				final float identity=identity(p);require(Float.isFinite(identity) && identity>=0 && identity<=1, "MODEL identity must be the finite observed primary or edge-retry alignment score");
				m.events[4]++;byPass[pass][4]++;if(identity>=cut){m.events[5]++;byPass[pass][5]++;}
			}else{
				require(p.termEquals("WINDOW", 0), "Unknown stage event type");m.events[2]++;byPass[pass][2]++;
				final int n=p.parseInt(10);require(n>=-1, "Unreached shortlist is -1; reached sizes are nonnegative");
				if(n>0){m.events[3]+=n;byPass[pass][3]+=n;}
				require(!p.termEquals("ACCEPT_HBM", 8) && !p.termEquals("ACCEPT_RESCUE", 8) && !p.termEquals("REJECT_RESCUE", 8), "A rescue outcome contradicts the declared raw equal-threshold arm");
				if(p.termStartsWith("ACCEPT", 8)){m.events[6]++;byPass[pass][6]++;}
			}
		}}finally{require(!in.close(), "Stage diagnostics read failed");}
	}
	private void calls(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			if(row.length==0 || row[0]=='#'){continue;}p.set(row);require(p.terms()==9, "Caller output must retain nine-column GFF");if(!p.termEquals("RNA", 2)){continue;}
			final String attrs=p.parseString(8);final int at=attrs.indexOf("model:");require(at>=0, "Generic RNA must retain an attributable model");
			final int start=at+6, end=start+prefix.length();
			if(!attrs.startsWith(prefix, start) || end<attrs.length() && attrs.charAt(end)!='_' && !Character.isWhitespace(attrs.charAt(end)) && attrs.charAt(end)!=',' && attrs.charAt(end)!=';'){continue;}
			final Member m=byId.get(RrnaResourceIO.first(p.parseString(0)));require(m!=null, "Emitted call outside the declared input population");
			require(p.parseInt(3)>0 && p.parseInt(4)>=p.parseInt(3) && p.parseInt(4)<=m.length && (p.termEquals("+", 6) || p.termEquals("-", 6)), "GFF bounds/strand must fit the input core");m.events[9]++;
		}}finally{require(!in.close(), "Final GFF read failed");}
	}
	private void write(String out){
		RrnaResourceIO.fresh(out, ".loci.tsv", ".summary.tsv", ".passes.tsv", ".called.ids.txt");
		final ByteStreamWriter loci=RrnaResourceIO.writer(out+".loci.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv"), passes=RrnaResourceIO.writer(out+".passes.tsv"), called=RrnaResourceIO.writer(out+".called.ids.txt");
		final long[] events=new long[NAMES.length], reached=new long[NAMES.length];final ByteBuilder b=new ByteBuilder();
		try{
			b.append("id\tlabel\tlength\tseed_observation");for(String name:NAMES){b.tab().append(name);}loci.print(b.append("\tHBM_EVALUATED\tHBM_PASSED\n"));
			for(Member m:members){
				b.clear().append(m.id).tab().append(label).tab().append(m.length).tab().append(m.events[0]==0?"NOT_OBSERVED":"OBSERVED");
				for(int i=0; i<NAMES.length; i++){events[i]+=m.events[i];if(m.events[i]>0){reached[i]++;}b.tab().append(m.events[i]);}
				loci.print(b.append("\tNA_NOT_OBSERVED\tNA_NOT_OBSERVED\n"));if(m.events[9]>0){called.println(m.id);}
			}
			summary.println("label\tfamily\tid_pass\tstage\tinput_records\tinput_bases\treached_records\tevents\treached_per_mb\tevents_per_mb\tinterpretation");
			for(int i=0; i<NAMES.length; i++){
				summary.print(b.clear().append(label).tab().append(family).tab().appendSlow(cut).tab().append(NAMES[i]).tab().append(members.size()).tab().append(bases)
					.tab().append(reached[i]).tab().append(events[i]).tab().append(reached[i]*1000000.0/bases, 8).tab().append(events[i]*1000000.0/bases, 8).append("\tRAW_CALLER_STAGE_ACCOUNTING_NOT_BIOLOGICAL_TRUTH\n"));
			}
			for(String name:new String[]{"HBM_EVALUATED", "HBM_PASSED"}){
				summary.print(b.clear().append(label).tab().append(family).tab().appendSlow(cut).tab().append(name).tab().append(members.size()).tab().append(bases)
					.append("\tNA\tNA\tNA\tNA\tNOT_OBSERVED_IN_STAGE_DIAGNOSTICS\n"));
			}
			passes.println("label\tpass\tstage\tevents");for(int pass=1; pass<=2; pass++){for(int i=2; i<=6; i++){
				passes.print(b.clear().append(label).tab().append(pass).tab().append(NAMES[i]).tab().append(byPass[pass][i]).nl());
			}}
		}finally{RrnaResourceIO.finish(loci, summary, passes, called);}
		System.err.println("NCRNA_SHRED_STAGES_PASS records="+members.size()+" emittedEvents="+events[9]+" emittedRecords="+reached[9]+"; HBM evaluation counts unobserved in these diagnostics");
	}
	/** appendSlow(float) may use exponents; retain zero-allocation parsing otherwise. */
	private static float identity(LineParser1 p){
		p.setBounds(13);final byte[] row=p.line();
		for(int i=p.a(); i<p.b(); i++){if(row[i]=='e' || row[i]=='E'){return Float.parseFloat(p.parseString(13));}}
		return p.parseFloat(13);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {Member(String i, int n){id=i;length=n;}final String id;final int length;final long[] events=new long[NAMES.length];}
	private final ArrayList<Member> members=new ArrayList<Member>();private final ObjectMap<String,Member> byId=new ObjectMap<String,Member>(String.class, Member.class);
	private final long[][] byPass=new long[3][NAMES.length];private long bases;private String family, prefix, label;private float cut=Float.NaN, border=Float.NaN;
	private static final String[] NAMES={"SEED_OBSERVED", "SEED_HIT", "WINDOW", "TERMINAL_SHORTLIST_ENTRY", "MODEL_ALIGNMENT_OBSERVED", "MODEL_IDENTITY_PASS_OBSERVED", "VERIFIER_ACCEPT", "DP_INPUT", "DP_RETAINED", "EMITTED"};
}
