package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

import prok.GradeNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaBenchmark;

/** Small saved-row oracles; invoked by the existing benchmark test entrypoint.
 * Equal subtype TP totals hide three losses and three gains. One lost parent
 * stays union-TP via another model, one is rematched within the same subtype,
 * and one becomes union-FN. Repeated parent/query IDs in a different input must
 * not collide. No calling, regrading or inferred reference labels.
 * @author Ganyu
 */
public final class NcrnaSubtypeComparisonTest{

	static void run(final File root){
		try{Files.createDirectory(root.toPath());}
		catch(IOException failure){throw new IllegalStateException("Cannot create fresh subtype fixture", failure);}
		final Fixture f=new Fixture(root);
		f.write();
		SummarizeNcrnaBenchmark.main(f.args("plain", false));
		SummarizeNcrnaBenchmark.main(f.args("joined", true));
		for(final String suffix:new String[]{"cells", "timing", "grades", "comparisons", "changes"}){
			require(lines(f.path("plain."+suffix+".tsv")).equals(lines(f.path("joined."+suffix+".tsv"))),
				"Optional subtype join changed ordinary output: "+suffix);
		}
		for(final String suffix:new String[]{"subtype_comparisons", "subtype_changes"}){
			require(!new File(f.path("plain."+suffix+".tsv")).exists(), "Absent optional inputs must not create subtype reports");
		}
		final ArrayList<String> summary=lines(f.path("joined.subtype_comparisons.tsv"));
		final ArrayList<String> changes=lines(f.path("joined.subtype_changes.tsv"));
		for(final String view:VIEWS){
			final String context="POLICY_VS_BASELINE\tpublic\t1\t1\tcandidate\t1\t1\t"+view+"\tRF01482\t";
			require(summary.contains(context+"4\t4\t3\t3\t1\t1\t1\t0\t1\t0\t2\t0\tSAVED_MATCH_WITNESSES_NONADDITIVE_NOT_BIOLOGICAL_SUBTYPE_PROOF"),
				"Equal subtype totals must expose three losses/gains and the independent union classifications");
			row(changes, context+"LOST\t"+A+"\tp1\tc\t100\t119\t+\tTP\tFN\tTP\t2\t"+MODEL1+
				"\tTP\t2\t"+MODEL0+"\tright\tOTHER_SUBTYPE_MODEL");
			row(changes, context+"LOST\t"+A+"\tp2\tc\t200\t219\t+\tTP\tFN\tTP\t3\t"+MODEL1+
				"\tTP\t3\t"+MODEL1+"\tright\tSAME_SUBTYPE_REMATCH");
			row(changes, context+"LOST\t"+A+"\tp4\tc\t400\t419\t+\tTP\tFN\tTP\t4\t"+MODEL1+
				"\tFN\t.\t.\tright\tUNION_MISS");
			row(changes, context+"GAINED\t"+A+"\tp7\tc\t800\t819\t+\tFN\tTP\tTP\t6\t"+MODEL0+
				"\tTP\t6\t"+MODEL1+"\tleft\tOTHER_SUBTYPE_MODEL");
			for(final String value:changes){
				require(!value.startsWith(context+"LOST\tB\t") && !value.startsWith(context+"GAINED\tB\t"),
					"Same parent/query IDs in input B must remain an independent unchanged input");
			}
		}
		require(changes.size()==13, "Two views each have exactly six changed parents; masks/unchanged inputs add none");

		final ArrayList<String> missing=lines(f.path("subtype_raw.tsv"));
		for(int i=missing.size()-1; i>0; i--){if(missing.get(i).startsWith("candidate\t")){missing.remove(i);}}
		f.writeRows("missing.tsv", missing);
		reject(f.modified("missing", "subtyperaw", f.path("missing.tsv")), "disagree with native summary counts");

		final ArrayList<String> duplicate=lines(f.path("subtype_raw.tsv"));
		for(final String line:lines(f.path("subtype_raw.tsv"))){
			if(line.contains("\tREFERENCE\tp1\t")){duplicate.add(line); break;}
		}
		f.writeRows("duplicate.tsv", duplicate);
		reject(f.modified("duplicate", "subtyperaw", f.path("duplicate.tsv")), "reference parent repeats");

		final ArrayList<String> badLine=lines(f.path("subtype_raw.tsv"));
		for(int i=1; i<badLine.size(); i++){
			if(badLine.get(i).contains("\tQUERY\t.\t2\t")){
				badLine.set(i, badLine.get(i).replace("\tQUERY\t.\t2\t", "\tQUERY\t.\t999\t"));
				break;
			}
		}
		f.writeRows("bad-line.tsv", badLine);
		reject(f.modified("bad-line", "subtyperaw", f.path("bad-line.tsv")), "no original ordinary query witness");

		final ArrayList<String> badModel=lines(f.path("subtype_raw.tsv"));
		for(int i=1; i<badModel.size(); i++){
			if(badModel.get(i).contains("\tQUERY\t")){
				badModel.set(i, badModel.get(i).replace("_m", "_different"));
				break;
			}
		}
		f.writeRows("bad-model.tsv", badModel);
		reject(f.modified("bad-model", "subtyperaw", f.path("bad-model.tsv")), "model or coordinates differ");

		final ArrayList<String> badParent=lines(f.path("subtype_raw.tsv"));
		for(int i=1; i<badParent.size(); i++){
			if(badParent.get(i).contains("\tREFERENCE\tp1\t")){
				badParent.set(i, badParent.get(i).replace("\tREFERENCE\tp1\t", "\tREFERENCE\tunknown\t"));
				break;
			}
		}
		f.writeRows("bad-parent.tsv", badParent);
		reject(f.modified("bad-parent", "subtyperaw", f.path("bad-parent.tsv")), "parent is absent");

		final String[] unpaired=f.args("unpaired", true);
		reject(Arrays.copyOf(unpaired, unpaired.length-1), "must be supplied together");
		final String[] timing=Arrays.copyOf(f.args("timing-conflict", true), unpaired.length+1);
		timing[timing.length-1]="timingonly=t";
		reject(timing, "does not accept subtype reports");
	}

	private static final class Fixture{
		Fixture(final File root_){root=root_;}
		String path(final String name){return new File(root, name).toString();}
		String[] args(final String name, final boolean subtype){
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+path("times.tsv"),
				"gradeinputs="+path("inputs.tsv"), "grades="+path("grades.tsv"), "raw="+path("raw.tsv"),
				"out="+path(name), "policies=public,candidate", "baseline=public", "workers=1", "repeats=1", "t=1"));
			if(subtype){
				args.add("subtyperaw="+path("subtype_raw.tsv"));
				args.add("subtypegrades="+path("subtype_summary.tsv"));
			}
			return args.toArray(new String[args.size()]);
		}
		String[] modified(final String name, final String option, final String value){
			final String[] args=args(name, true);
			for(int i=0; i<args.length; i++){if(args[i].startsWith(option+"=")){args[i]=option+"="+value;}}
			return args;
		}
		void write(){
			try(Output times=new Output(path("times.tsv")); Output inputs=new Output(path("inputs.tsv"));
				Output grades=new Output(path("grades.tsv")); Output raw=new Output(path("raw.tsv"));
				Output sr=new Output(path("subtype_raw.tsv")); Output ss=new Output(path("subtype_summary.tsv"))){
				text(times, "policy\tthreads\trepeat\tuser_seconds\tsystem_seconds\twall_seconds");
				text(inputs, CONTEXT+GradeNcrnaReference.INPUT_HEADER);
				text(grades, CONTEXT+GradeNcrnaReference.SUMMARY_HEADER);
				text(raw, CONTEXT+GradeNcrnaReference.RAW_HEADER);
				text(sr, CONTEXT+SUBTYPE_RAW_HEADER);
				text(ss, CONTEXT+SUBTYPE_SUMMARY_HEADER);
				for(final String policy:new String[]{"public", "candidate"}){
					final boolean changed=policy.equals("candidate");
					final String context=policy+"\t1\t1\t";
					text(times, context+"12\t1\t14");
					text(inputs, context+A+"\t/"+A+".fa\t/qA\t/rA\t4\t0\t0\t0\t4\t0");
					text(inputs, context+"B\t/B.fa\t/qB\t/rB\t1\t0\t0\t0\t1\t0");
					text(inputs, context+"empty\t/empty.fa\t/qE\t/rE\t0\t0\t0\t0\t0\t0");
					text(grades, context+"cobalamin\t3\t5\t0\t2\t7\t1\t0\t0\t0\t5\t5\t5\t1\t0.714286\t0.5\tFIXTURE");
					for(int i=0; i<PARENTS.length; i++){
						final int line=(changed ? CANDIDATE_LINES : PUBLIC_LINES)[i];
						final String input=i==5 ? "B" : A;
						final String shape="c\t"+STARTS[i]+"\t"+(STARTS[i]+19)+"\t+\t";
						final String state=i==7 ? "DONT_CARE" : line>0 ? "TP" : "FN";
						text(raw, context+input+"\tREFERENCE\tcobalamin\t"+shape+state+"\t"+PARENTS[i]+
							"\t.\t.\t"+(line>0 ? line : ".")+"\t.\t.\t.\tRF01482\t.\t.");
						if(line>0){
							final String model=model(changed, i);
							text(raw, context+input+"\tQUERY\tcobalamin\t"+shape+"TP\t.\t"+line+"\t"+PARENTS[i]+
								"\t.\t.\tfalse\ttrue\t"+model+"\t0\t0");
						}
					}
					for(final String view:VIEWS){
						for(final String type:TYPES){
							// Literal independently counted totals: RF01482 has four TPs
							// in each build but three nonnested losses and three gains.
							final String totals=type.equals("RF01482") ? "3\t4\t0\t3\t7\t1\t0\t0\t"+
								(changed ? "3\t3\t3" : "4\t4\t4")+"\t1\t0.571429"
								: type.equals("RF00174") ? "3\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA"
								: "3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA";
							text(ss, context+view+"\t"+type+"\t"+totals+"\t0.5\tFIXTURE");
							if(type.equals("RF01482")){
								for(int i=0; i<PARENTS.length; i++){
									final int line=(changed ? CANDIDATE_SUBTYPE_LINES : PUBLIC_SUBTYPE_LINES)[i];
									final String input=i==5 ? "B" : A;
									final String state=i==7 ? "DONT_CARE" : line>0 ? "TP" : "FN";
									text(sr, context+view+"\t"+type+"\t"+input+"\tREFERENCE\t"+PARENTS[i]+
										"\t.\tRF01482\tc\t"+STARTS[i]+"\t"+(STARTS[i]+19)+"\t+\t"+state+
										"\t.\t"+(line>0 ? line : ".")+"\t.\t.\t.\tRF01482\t.\tfalse\t.\t.");
								}
							}
							for(int i=0; i<PARENTS.length; i++){
								final int line=(changed ? CANDIDATE_LINES : PUBLIC_LINES)[i];
								if(line<0){continue;}
								final String model=model(changed, i);
								if(!model.startsWith("cobalamin_"+type)){continue;}
								final String input=i==5 ? "B" : A;
								final boolean matched=type.equals("RF01482");
								final String parent=changed && i==1 ? "p3" : PARENTS[i];
								final String endpoints=!matched ? ".\t." : changed && i==1 ? "-1\t-1" : "0\t0";
								text(sr, context+view+"\t"+type+"\t"+input+"\tQUERY\t.\t"+line+"\t"+model+
									"\tc\t"+STARTS[i]+"\t"+(STARTS[i]+19)+"\t+\t"+(matched ? "TP" : "FP")+
									"\t"+(matched ? parent : ".")+"\t.\t.\tfalse\t"+matched+"\t"+(matched ? "RF01482" : ".")+
									"\t.\t"+(matched ? "false" : ".")+"\t"+endpoints);
							}
						}
					}
				}
			}
		}
		void writeRows(final String name, final ArrayList<String> rows){
			try(Output out=new Output(path(name))){for(final String row:rows){text(out, row);}}
		}
		final File root;
	}

	private static String model(final boolean candidate, final int parent){
		return (candidate && parent==0) || (!candidate && parent==6) ? MODEL0 : MODEL1;
	}
	private static void text(final Output out, final String row){out.writer.println(row.getBytes(StandardCharsets.UTF_8));}
	private static ArrayList<String> lines(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result;
	}
	private static void row(final ArrayList<String> rows, final String expected){
		require(rows.contains(expected), "Missing exact subtype witness row: "+expected);
	}
	private static void reject(final String[] args, final String diagnostic){
		try{SummarizeNcrnaBenchmark.main(args);}
		catch(IllegalArgumentException failure){
			require(failure.getMessage().contains(diagnostic), "Unexpected rejection: "+failure.getMessage());
			return;
		}
		throw new AssertionError("Invalid subtype input was accepted: "+diagnostic);
	}
	private static void require(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private static final String A="α%25";
	private static final String MODEL0="cobalamin_RF00174_m";
	private static final String MODEL1="cobalamin_RF01482_m";
	private static final String[] VIEWS={"MODEL_SUPPORT_RECOVERY", "EXCLUSIVE_SUBTYPE_ACCURACY"};
	private static final String[] TYPES={"RF00174", "RF01482", "RF01689"};
	private static final String[] PARENTS={"p1", "p2", "p3", "p4", "p5", "p1", "p7", "mask"};
	private static final int[] STARTS={100, 200, 201, 400, 500, 100, 800, 900};
	private static final int[] PUBLIC_LINES={2, 3, -1, 4, -1, 2, 6, -1};
	private static final int[] CANDIDATE_LINES={2, 3, -1, -1, 5, 2, 6, -1};
	private static final int[] PUBLIC_SUBTYPE_LINES={2, 3, -1, 4, -1, 2, -1, -1};
	private static final int[] CANDIDATE_SUBTYPE_LINES={-1, -1, 3, -1, 5, 2, 6, -1};
	private static final String CONTEXT="policy\tthreads\trepeat\t";
	private static final String SUBTYPE_RAW_HEADER="view\tsubtype\tinput_id\trole\tparent_id\tquery_line\tmodel\tseqid\tstart\tstop\tstrand\toutcome\tmatched_parent_id\tmatched_query_line\tmask_parent_id\tduplicate\tscored_eligible\tscored_support\tmarginal_support\tmulti_supported\td5\td3";
	private static final String SUBTYPE_SUMMARY_HEADER="view\tsubtype\tinputs\tTP\tFP\tFN\tscored_references\tmasked_references\tmasked_calls\tduplicate_queries\texact_5prime\texact_3prime\texact_both\tprecision\trecall\toverlap\tqualification";
}
