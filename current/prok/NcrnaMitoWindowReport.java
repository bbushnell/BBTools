package prok;

import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Cap-two work projection on observed windows, never a new caller recall estimate.
 * @author Brian, Ganyu */
public final class NcrnaMitoWindowReport {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaMitoWindowReport.class, false);
		try{final NcrnaMitoWindowReport x=new NcrnaMitoWindowReport();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] pair=Parse.splitOnFirst(arg, '=');require(pair.length==2, "Use flag=value");
			final String key=pair[0].toLowerCase(Locale.ROOT);require(KNOWN.contains("|"+key+"|") && flags.put(key, pair[1])==null, "Unknown/repeated flag: "+key);}
		in=required("in");time=required("time");prefix=required("outprefix");label=required("label");require(label.matches("[A-Za-z0-9_.-]+"), "Label must be a safe table token");
		final ObjectMap<String, String> pins=EukPanelGffValidation.hashes(required("verified"));hash=pins.get(in);
		require(hash!=null && SeedOffsetTable.validSha80(hash), "Window input must have a producer-verified sha80 binding");
		timeHash=pins.get(time);require(timeHash!=null && SeedOffsetTable.validSha80(timeHash), "Timing input must have a producer-verified sha80 binding");
		RrnaResourceIO.fresh(prefix, ".summary.tsv", ".inputs.tsv", ".top_shared.tsv");Shared.setThreads(1);ReadWrite.setZipThreads(1);
	}
	private void run(){
		readTime();
		final ByteFile reader=RrnaResourceIO.open(in);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(reader, MitoSsuTierLog.HEADER, in);for(byte[] line;(line=reader.nextLine())!=null;){
			p.set(line);require(p.terms()==18, "Observer rows require all eighteen fields");
			require(p.setBounds(0)>0, "Every observed window needs a contig name");
			final long strand=number(p, 1), pass=number(p, 2), start=number(p, 3), stop=number(p, 4), unique=number(p, 5);
			final long top=number(p, 6), second=number(p, 7), third=number(p, 8), models=number(p, 9), calls=number(p, 10), skips=number(p, 11);
			final boolean accepted=p.termEquals("ACCEPT", 16);
			final long model=signedId(p, 12), rank=signedId(p, 13), attempt=number(p, 14), shared=number(p, 15), firstTwo=number(p, 17);
			require(strand<=1 && pass>0 && start<=stop && stop<=Integer.MAX_VALUE && unique<=Integer.MAX_VALUE, "Window coordinates, strand and query count must remain valid");
			require(top>=second && second>=third && top<=65535, "Shared counts must retain native descending K9 ranking");
			allTop[(int)top]=add(allTop[(int)top], 1);
			require(models<=Integer.MAX_VALUE && calls>=models && calls<=2*models, "Each attempted model has one full verifier and at most one physical-edge retry");
			require(firstTwo>=Math.min(2, models) && firstTwo<=2*Math.min(2, models) && firstTwo<=calls
				&& (models>2 || firstTwo==calls), "First-two work must be the actual bounded prefix of this window's verifier calls");
			require(accepted || p.termEquals("REJECT_VERIFY", 16) || p.termEquals("REJECT_ENTRY", 16)
				|| p.termEquals("REJECT_RELATIVE", 16) || p.termEquals("REJECT_CAP", 16), "Unknown observer outcome");
			if(accepted){
				require(model>=0 && rank>=attempt && attempt>0 && attempt==models && shared>0 && shared<=top, "An observed acceptance must bind the last attempted model and its ranking");
				accepts=add(accepts, 1);minimumShared=Math.min(minimumShared, shared);acceptedTop[(int)top]=add(acceptedTop[(int)top], 1);
				if(attempt<=2){retained=add(retained, 1);minimumRetainedShared=Math.min(minimumRetainedShared, shared);}
				else{late=add(late, 1);}
			}else{require(model==-1 && rank==-1 && attempt==0 && shared==0, "Rejected windows cannot carry invented acceptance metadata");}
			windows=add(windows, 1);attempts=add(attempts, models);alignments=add(alignments, calls);lengthSkips=add(lengthSkips, skips);capCalls=add(capCalls, firstTwo);
			if(models>2){beyondTwo=add(beyondTwo, 1);}if(models==0){zeroAttempt=add(zeroAttempt, 1);}
		}}finally{require(!reader.close(), "Observer input read failed");}
		require(accepts==retained+late && capCalls<=alignments, "Cap-two projection must conserve observed acceptances and cannot add work");
		write();
	}
	private void readTime(){
		final ByteFile reader=RrnaResourceIO.open(time);final LineParser1 p=new LineParser1(' ');
		try{final byte[] line=reader.nextLine();require(line!=null, "Missing run timing");p.set(line);require(p.terms()==3 && reader.nextLine()==null, "Require one wall/user/system timing row");
			final String[] keys={"wall_seconds=", "user_seconds=", "system_seconds="};
			for(int i=0;i<3;i++){final String term=p.parseString(i);require(term.startsWith(keys[i]), "Unexpected timing field");seconds[i]=Double.parseDouble(term.substring(keys[i].length()));
				require(Double.isFinite(seconds[i]) && seconds[i]>=0, "Run timing must be finite and nonnegative");}
			require(Double.isFinite(seconds[1]+seconds[2]), "CPU seconds cannot overflow");
		}finally{require(!reader.close(), "Timing input read failed");}
	}
	private void write(){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".summary.tsv"), inputs=RrnaResourceIO.writer(prefix+".inputs.tsv"), hist=RrnaResourceIO.writer(prefix+".top_shared.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println(HEADER);b.append(label).tab().append(windows).tab().append(accepts).tab().append(windows-accepts)
			.tab().append(retained).tab().append(late).tab().append(attempts).tab().append(alignments).tab().append(capCalls)
			.tab().append(alignments-capCalls).tab().append(lengthSkips).tab().append(beyondTwo).tab().append(zeroAttempt).tab();
			if(accepts==0){b.append("NA");}else{b.append(minimumShared);}b.tab();if(retained==0){b.append("NA");}else{b.append(minimumRetainedShared);}
			for(double value:seconds){b.tab().append(value, 2);}b.tab().append(seconds[1]+seconds[2], 2);
			out.print(b.append("\tOBSERVED_WINDOWS_ONLY_NOT_CHANGED_RUN_RECALL_OR_TOTAL_WORK\n"));
			inputs.println("role\tpath\tsha80");inputs.print(b.clear().append("windows\t").append(in).tab().append(hash).nl());
			inputs.print(b.clear().append("timing\t").append(time).tab().append(timeHash).nl());
			hist.println("label\ttop_shared9\tall_entered_window_events\taccepted_window_events\tinterpretation");
			for(int score=0;score<allTop.length;score++){if(allTop[score]>0){hist.print(b.clear().append(label).tab().append(score).tab().append(allTop[score]).tab().append(acceptedTop[score]).append("\tWINDOW_EVENTS_NOT_UNIQUE_LOCI\n"));}}
		}finally{RrnaResourceIO.finish(out, inputs, hist);}
		System.err.println("MITO_WINDOW_REPORT_PASS windows="+windows+" accepted="+accepts+" cap2_retained="+retained+" new_alignments=0");
	}
	private static long number(LineParser1 p, int term){
		p.setBounds(term);require(p.a()<p.b(), "Missing numeric observer field");final byte[] line=p.line();long value=0;
		for(int i=p.a();i<p.b();i++){final int digit=line[i]-'0';require(digit>=0 && digit<=9 && value<=(Long.MAX_VALUE-digit)/10, "Observer counts must be exact nonnegative long integers");value=value*10+digit;}return value;
	}
	private static long signedId(LineParser1 p, int term){return p.termEquals("-1", term)?-1:number(p, term);}
	private static long add(long a, long b){require(b>=0 && a<=Long.MAX_VALUE-b, "Observer totals must not overflow");return a+b;}
	private String required(String key){final String v=flags.get(key);require(v!=null && !v.isEmpty(), "Missing flag: "+key);return v;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String HEADER="label\twindows\tobserved_accepts\tobserved_rejects\tcap2_retained_accepts\tcap2_lost_late_accepts\tmodels_attempted\taligner_calls\tcap2_prefix_aligner_calls\tcap2_removed_aligner_calls\tlength_skips\twindows_beyond_two_models\twindows_zero_attempts\tminimum_accepted_shared9\tminimum_retained_shared9\twall_seconds\tuser_seconds\tsystem_seconds\tcpu_seconds\tinterpretation";
	private static final String KNOWN="|in|time|verified|outprefix|label|";
	private final ObjectMap<String, String> flags=new ObjectMap<String, String>(String.class, String.class);
	private String in, time, prefix, label, hash, timeHash;
	private final double[] seconds=new double[3];
	private long windows, accepts, retained, late, attempts, alignments, capCalls, lengthSkips, beyondTwo, zeroAttempt;
	private long minimumShared=Long.MAX_VALUE, minimumRetainedShared=Long.MAX_VALUE;
	private final long[] allTop=new long[65536], acceptedTop=new long[65536];
}
