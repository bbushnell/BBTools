package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectSet;
import map.IntLongHashMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Learns integer per-model endpoint medians and replays saved matched calls.
 * Adapted from template.A_SampleByteFile, using native stream owners for multiple
 * inputs. One thread owns every parser, histogram and counter. No sequences are
 * aligned, calls rematched, networks fitted, or detections measured here.
 * Training excludes exact assemblies from the supplied outer-panel audit.
 * Quantiles use nearest rank (the lower median for even n), as NcrnaEndpointSummary.
 * Replay holds both ends if shifted one-based coordinates become invalid. Contig
 * lengths and gene-path competition are unavailable and are not reconstructed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelEndpointOffsets{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaModelEndpointOffsets.class, false);
		try{
			new NcrnaModelEndpointOffsets(pp.args).process();
			pp.outstream.println("NCRNA_MODEL_ENDPOINT_OFFSETS_COMPLETE_NO_CALLER_NO_FIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Requires complete saved populations, an explicit catalog and fresh outputs. */
	public NcrnaModelEndpointOffsets(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("refs")){refs=b.split(",", -1);}
			else if(a.equals("compare")){comparison=b.split(",", -1);}
			else if(a.equals("members")){members=b.split(",", -1);}
			else if(a.equals("family")){family=b;}
			else if(a.equals("evaluation")){evaluation=b.split(",", -1);}
			else if(a.equals("panels")){panels=b.split(",", -1);}
			else if(a.equals("bindings")){bindings=b;}
			else if(a.equals("exclude")){exclude=b;}
			else if(a.equals("expectedmodels")){expectedModels=Integer.parseInt(b);}
			else if(a.equals("expected")){expected=counts(b);}
			else if(a.equals("breakdown")){breakdownMember=b;}
			else if(a.equals("breakdowncounts")){breakdownCounts=counts(b);}
			else if(a.equals("centeredhistograms")){centeredHistograms=Parse.parseBoolean(b);}
			else if(a.equals("refit")){refit=Parse.parseBoolean(b);}
			else if(a.equals("minlengths")){minimumLengths=counts(b);}
			else if(a.equals("thresholds")){
				final long[] values=counts(b);
				thresholds=new int[values.length];
				for(int i=0; i<values.length; i++){
					require(values[i]>0 && values[i]<=Integer.MAX_VALUE && (i==0 || values[i]>values[i-1]),
						"Support thresholds must be distinct, positive and ascending");
					thresholds[i]=(int)values[i];
				}
			}else{require(parser.parse(arg, a, b), "Unknown endpoint-offset argument: "+arg);}
		}
		training=parser.in1==null ? null : parser.in1.split(",", -1);
		out=parser.out1;
		if(comparison!=null){
			require(comparison.length==2 && out!=null && members!=null && members.length>0
				&& family!=null && !family.isEmpty() && expectedModels>0 && thresholds.length==1,
				"compare=old.models.tsv,new.models.tsv requires members= family= expectedmodels= thresholds=one_minimum out=prefix");
			require(training==null && refs==null && evaluation==null && panels==null && bindings==null
				&& exclude==null && expected==null && breakdownMember==null && breakdownCounts==null && !centeredHistograms && minimumLengths==null && !refit,
				"Saved-table comparison never reads witnesses, references or evaluation panels");
			require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
				&& parser.maxReads<0 && parser.samplerate==1, "Compare complete saved tables without overwrite or sampling");
			final ObjectSet<String> labels=new ObjectSet<String>(String.class);
			labels.add(family);
			for(final String member:members){require(!member.isEmpty() && labels.add(member), "Unique member and union labels are required");}
			require(Tools.testInputFiles(false, true, comparison), "Both saved model tables must be readable and distinct");
			require(Tools.testOutputFiles(false, false, false, out+".stability.tsv", out+".stability_summary.tsv"), "Stability outputs must be fresh");
			return;
		}
		if(refit){
			require(centeredHistograms && thresholds.length==1 && exclude==null && evaluation==null && panels==null
				&& breakdownMember==null && breakdownCounts==null && expected!=null && expected.length==2 && expected[0]==expected[1],
				"All-data refit requires centeredhistograms=t, one support threshold, equal selected/retained counts and no exclusions or evaluation");
			evaluation=new String[0];
			panels=new String[0];
		}
		require(training!=null && out!=null && refs!=null && members!=null && family!=null
			&& evaluation!=null && panels!=null && bindings!=null && (refit || exclude!=null) && expected!=null && expectedModels>0,
			"Require in=training_witnesses refs= members= family= evaluation= panels= bindings= exclude= expected= expectedmodels= out=prefix");
		require(refs.length==members.length && members.length>0 && evaluation.length==panels.length
			&& (refit || evaluation.length>0) && expected.length==evaluation.length+2,
			"Expected counts are training-before-exclusion, training-retained, then one per evaluation panel");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Replay complete populations without overwrite, append or sampling");
		final ObjectSet<String> labels=new ObjectSet<String>(String.class);
		labels.add(family);
		for(final String member:members){require(!member.isEmpty() && labels.add(member), "Member and union labels must be unique");}
		labels.clear();
		for(final String panel:panels){require(!panel.isEmpty() && labels.add(panel), "Panel labels must be unique and nonempty");}
		if(breakdownMember!=null){
			for(int i=0; i<members.length; i++){if(members[i].equals(breakdownMember)){breakdownIndex=i;}}
			require(breakdownIndex>=0 && thresholds.length==1 && breakdownCounts!=null && breakdownCounts.length==panels.length,
				"breakdown=member requires one explicit support threshold and breakdowncounts=one accepted member count per panel");
			require(Tools.testOutputFiles(false, false, false, out+".breakdown.tsv"), "Breakdown output must be fresh and writable");
		}else{require(breakdownCounts==null, "breakdowncounts= requires breakdown=member");}
		if(minimumLengths!=null){
			require(centeredHistograms && minimumLengths.length==members.length,
				"minlengths= requires centeredhistograms=t and one explicit family minimum per ordered member");
			for(final long minimum:minimumLengths){
				require(minimum>0 && minimum<=Integer.MAX_VALUE, "Each bound family minimum must be a positive integer length");
			}
		}
		if(centeredHistograms){
			require(breakdownMember==null && thresholds.length==1,
				"Centered histogram export requires one support threshold and cannot also request a panel breakdown");
			for(final String member:members){
				require(member.matches("[A-Za-z0-9_]+"), "Histogram member labels must be safe literal filename components");
				require(Tools.testOutputFiles(false, false, false, out+"."+member+".histogram.tsv"), "Centered histograms must be fresh");
			}
			require(Tools.testOutputFiles(false, false, false, out+".offsets.tsv", out+".training.tsv"), "Offset resource and training accounting must be fresh");
		}
		for(final String path:training){require(Tools.testInputFiles(false, true, path), "Unreadable training witnesses: "+path);}
		for(final String path:evaluation){require(Tools.testInputFiles(false, true, path), "Unreadable panel witnesses: "+path);}
		require(Tools.testInputFiles(false, true, refs) && Tools.testInputFiles(false, true, bindings, exclude), "Unreadable catalog or source exclusions");
		require(Tools.testOutputFiles(false, false, false, out+".models.tsv", out+".support.tsv", out+".single.tsv", out+".both.tsv", out+".population.tsv"),
			"All replay output files must be fresh and writable");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads all inputs and checks accepted counts before publishing any result. */
	public void process(){
		if(comparison!=null){
			compareModels();
			return;
		}
		loadModels();
		final Sources sources=new Sources(bindings, exclude);
		long selected=0;
		long retained=0;
		final ObjectSet<Target> seenReferences=new ObjectSet<Target>(Target.class);
		final ObjectSet<Target> seenCalls=new ObjectSet<Target>(Target.class);
		for(final String path:training){
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, HEADER, path);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					readRow(line, seenReferences, seenCalls);
					selected++;
					if(sources.excludes(row)){continue;}
					retained++;
					models.get(modelIndex).errors[0].add(error5);
					models.get(modelIndex).errors[1].add(error3);
				}
			}
		}
		require(selected==expected[0] && retained==expected[1], "Training witness counts disagree with accepted producer counts: "+selected+"/"+retained);
		for(final Model model:models){model.finish();}
		if(centeredHistograms){
			writeCenteredHistograms(sources);
			if(refit){
				writeModels();
				writeSupport();
				try(Output population=new Output(out+".population.tsv")){
					population.writer.println("training_selected\ttraining_retained\ttraining_excluded\tmodels\tquantiles\tqualification");
					population.writer.print(buffer.clear().append(selected).tab().append(retained).append("\t0\t")
						.append(models.size()).append("\tNEAREST_RANK_LOWER_MEDIAN\tALL_DATA_REFIT_NO_HELDOUT_NO_EVALUATION\n"));
				}
			}
			return;
		}
		if(breakdownMember!=null){writeBreakdown();return;}
		final Stats[][][][] results=new Stats[evaluation.length][thresholds.length][members.length+1][3];
		for(int panel=0; panel<evaluation.length; panel++){replay(panel, results[panel]);}
		writeModels();
		writeSupport();
		try(Output single=new Output(out+".single.tsv"); Output both=new Output(out+".both.tsv");
				Output population=new Output(out+".population.tsv")){
			single.writer.println(RESULT_HEADER);
			both.writer.println(RESULT_HEADER);
			population.writer.println("training_selected\ttraining_retained\ttraining_excluded\tmodels\tquantiles\tqualification");
			population.writer.print(buffer.clear().append(selected).tab().append(retained).tab().append(selected-retained)
				.tab().append(models.size()).append("\tNEAREST_RANK_LOWER_MEDIAN\t").append(QUALIFICATION).nl());
			for(int panel=0; panel<evaluation.length; panel++){
				for(int t=0; t<thresholds.length; t++){
					for(int member=0; member<=members.length; member++){
						for(int end=0; end<3; end++){
							final Stats stats=results[panel][t][member][end];
							stats.row(buffer.clear(), panels[panel], label(member), thresholds[t], end);
							(end==2 ? both : single).writer.print(buffer);
						}
					}
				}
			}
		}
	}

	/** Keeps every catalog model, including those with no training match. */
	private void loadModels(){
		for(int member=0; member<refs.length; member++){
			final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(refs[member]);
			for(final String name:library.names){
				final byte[] bytes=name.getBytes(java.nio.charset.StandardCharsets.UTF_8);
				probe.set(bytes, 0, bytes.length);
				require(!byModel.contains(probe), "Model identifiers must be unique across the full catalog: "+name);
				byModel.put(probe.copy(), models.size());
				models.add(new Model(name, member));
			}
		}
		require(models.size()==expectedModels, "Full model catalog count differs: expected="+expectedModels+" actual="+models.size());
	}

	/** Validates actual saved winner and transcript-oriented signed errors. */
	private void readRow(final byte[] line, final ObjectSet<Target> seenReferences, final ObjectSet<Target> seenCalls){
		row.set(line);
		require(row.terms()==12, "Witness rows require all twelve exporter fields");
		for(int i=0; i<12; i++){require(row.length(i)>0, "Witness fields cannot be empty");}
		row.setBounds(3);
		int end=row.a();
		while(end<row.b() && (line[end]&255)>32){end++;}
		probe.set(line, row.a(), end);
		modelIndex=byModel.get(probe);
		require(modelIndex>=0, "Saved winning model is absent from the frozen catalog");
		require(row.termEquals(members[models.get(modelIndex).member], 2), "Saved winner and emitted member disagree");
		callStart=row.parseInt(5);
		callStop=row.parseInt(6);
		final int truthStart=row.parseInt(7);
		final int truthStop=row.parseInt(8);
		minus=row.termEquals("-", 9);
		require((minus || row.termEquals("+", 9)) && callStart>=1 && callStop>=callStart
			&& truthStart>=1 && truthStop>=truthStart, "Witnesses need valid one-based inclusive coordinates and strand");
		error5=minus ? truthStop-callStop : callStart-truthStart;
		error3=minus ? truthStart-callStart : callStop-truthStop;
		require(row.parseInt(10)==error5 && row.parseInt(11)==error3, "Signed errors must equal caller minus reference in transcript orientation");
		final long overlap=Math.max(0L, (long)Math.min(callStop, truthStop)-Math.max(callStart, truthStart)+1);
		require(2*overlap>=Math.max((long)callStop-callStart+1, (long)truthStop-truthStart+1), "Saved matched pair violates reciprocal-half overlap");
		buffer.clear();
		field(buffer, row, 0).tab();
		field(buffer, row, 1);
		probe.set(buffer.array, 0, buffer.length());
		require(seenReferences.add(probe.copy()), "A reference ID cannot contribute multiple matched calls");
		buffer.clear();
		field(buffer, row, 0).tab();
		field(buffer, row, 4).tab().append(callStart).tab().append(callStop).tab().append(minus ? '-' : '+');
		probe.set(buffer.array, 0, buffer.length());
		require(seenCalls.add(probe.copy()), "One physical call cannot consume multiple references");
	}

	/** Applies training-only offsets to fixed panel matches; never changes the cohort. */
	private void replay(final int panel, final Stats[][][] result){
		for(int t=0; t<thresholds.length; t++){
			for(int member=0; member<=members.length; member++){
				for(int end=0; end<3; end++){result[t][member][end]=new Stats();}
			}
		}
		final ObjectSet<Target> seenReferences=new ObjectSet<Target>(Target.class);
		final ObjectSet<Target> seenCalls=new ObjectSet<Target>(Target.class);
		long count=0;
		try(Input input=new Input(evaluation[panel])){
			RrnaResourceIO.header(input.reader, HEADER, evaluation[panel]);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				readRow(line, seenReferences, seenCalls);
				count++;
				final Model model=models.get(modelIndex);
				for(int t=0; t<thresholds.length; t++){
					final boolean supported=model.errors[0].size>=thresholds[t];
					final int median5=supported ? model.median[0] : 0;
					final int median3=supported ? model.median[1] : 0;
					final long shiftedStart=(long)callStart+(minus ? median3 : -(long)median5);
					final long shiftedStop=(long)callStop+(minus ? median5 : -(long)median3);
					final boolean invalid=shiftedStart<1 || shiftedStop<shiftedStart || shiftedStop>Integer.MAX_VALUE;
					final long after5=invalid ? error5 : (long)error5-median5;
					final long after3=invalid ? error3 : (long)error3-median3;
					for(int k=0; k<2; k++){
						final Stats[] stats=result[t][k==0 ? model.member : members.length];
						stats[0].accept(error5, 0, after5, 0, supported, invalid);
						stats[1].accept(error3, 0, after3, 0, supported, invalid);
						stats[2].accept(error5, error3, after5, after3, supported, invalid);
					}
				}
			}
		}
		require(count==expected[panel+2], "Panel witness count differs for "+panels[panel]+": expected="+expected[panel+2]+" actual="+count);
	}

	/** Compares nested training populations on the identical ordered model axis.
	 * Missing medians remain NA. Stability bins include only models supported in
	 * BOTH populations; newly supported models are counted separately.
	 */
	private void compareModels(){
		final Stability[][] totals=new Stability[members.length+1][2];
		for(int i=0; i<totals.length; i++){
			for(int end=0; end<2; end++){totals[i][end]=new Stability();}
		}
		final LineParser1 before=new LineParser1('\t');
		final LineParser1 after=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final ByteBuilder rows=new ByteBuilder("member\tmodel\tend\tn_before\tn_after\tmedian_before_bp\tmedian_after_bp\tmedian_change_bp\tiqr_before_bp\tiqr_after_bp\tsupported_before\tsupported_after\n");
		String lastName=null;
		int lastMember=-1;
		long lastBefore=-1;
		long lastAfter=-1;
		int count=0;
		try(Input oldTable=new Input(comparison[0]); Input newTable=new Input(comparison[1])){
			RrnaResourceIO.header(oldTable.reader, MODEL_HEADER, comparison[0]);
			RrnaResourceIO.header(newTable.reader, MODEL_HEADER, comparison[1]);
			for(byte[] line=oldTable.reader.nextLine(); line!=null; line=oldTable.reader.nextLine()){
				require(count<2L*expectedModels, "Saved tables exceed the declared complete model axis");
				final byte[] next=newTable.reader.nextLine();
				require(next!=null, "New model table ended before the original model axis");
				before.set(line);
				after.set(next);
				require(before.terms()==10 && after.terms()==10, "Saved model tables require the complete ten-column schema");
				for(int i=0; i<3; i++){
					before.setBounds(i);
					after.setBounds(i);
					require(before.length(i)>0 && before.length(i)==after.length(i), "Model identity widths differ between support populations");
					for(int j=0; j<before.length(i); j++){
						require(line[before.a()+j]==next[after.a()+j], "Member, model and end must match in their original order");
					}
				}
				final int end=count&1;
				require(before.termEquals(end==0 ? "5prime" : "3prime", 2), "Each catalog model requires ordered 5prime and 3prime rows");
				int member=-1;
				for(int i=0; i<members.length; i++){if(before.termEquals(members[i], 0)){member=i;}}
				require(member>=0, "Saved model belongs to an undeclared member");
				final long n0=before.parseLong(3);
				final long n1=after.parseLong(3);
				require(n0>=0 && n1>=n0, "Pooled support includes every old retained match; counts cannot decrease");
				if(end==0){
					lastName=before.parseString(1);
					require(names.add(lastName), "Every catalog model must appear exactly once per end");
					lastMember=member;
					lastBefore=n0;
					lastAfter=n1;
				}else{
					require(before.termEquals(lastName, 1) && member==lastMember && n0==lastBefore && n1==lastAfter,
						"Both ends of a model must share its name, member and matched-call support");
				}
				final long median0=savedMedian(before, n0);
				final long median1=savedMedian(after, n1);
				final boolean supported0=n0>=thresholds[0];
				final boolean supported1=n1>=thresholds[0];
				for(int i=0; i<3; i++){field(rows, before, i).tab();}
				rows.append(n0).tab().append(n1).tab();
				field(rows, before, 4).tab();
				field(rows, after, 4).tab();
				if(n0==0 || n1==0){rows.append("NA");}else{rows.append(median1-median0);}
				field(rows.tab(), before, 7).tab();
				field(rows, after, 7).tab().append(supported0 ? 1 : 0).tab().append(supported1 ? 1 : 0).nl();
				totals[member][end].add(supported0, supported1, median1-median0);
				totals[members.length][end].add(supported0, supported1, median1-median0);
				count++;
			}
			require(newTable.reader.nextLine()==null && count==2L*expectedModels, "Both saved tables must cover the complete declared model axis");
		}
		try(Output detail=new Output(out+".stability.tsv"); Output summary=new Output(out+".stability_summary.tsv")){
			detail.writer.print(rows);
			summary.writer.println("cohort\tend\tminimum_n\tmodels\tsupported_before\tsupported_after\tsupported_both\tnewly_supported\tmedian_unchanged\tmedian_change_1_2bp\tmedian_change_3_5bp\tmedian_change_above5bp");
			for(int member=0; member<totals.length; member++){
				for(int end=0; end<2; end++){
					final Stability s=totals[member][end];
					summary.writer.print(buffer.clear().append(label(member)).tab().append(end==0 ? "5prime" : "3prime")
						.tab().append(thresholds[0]).tab().append(s.models).tab().append(s.before).tab().append(s.after)
						.tab().append(s.both).tab().append(s.gained).tab().append(s.same).tab().append(s.small)
						.tab().append(s.moderate).tab().append(s.large).nl());
				}
			}
		}
	}

	/** Validates missing support and the ordered signed quantiles before comparison. */
	private static long savedMedian(final LineParser1 p, final long n){
		if(n==0){
			for(int i=4; i<=7; i++){require(p.termEquals("NA", i), "Zero support has no learned median or spread");}
			return 0;
		}
		final long median=p.parseInt(4);
		final long q25=p.parseInt(5);
		final long q75=p.parseInt(6);
		require(q25<=median && median<=q75 && p.parseLong(7)==q75-q25, "Signed quartiles must bracket the median and conserve IQR");
		return median;
	}

	private void writeModels(){
		try(Output output=new Output(out+".models.tsv")){
			output.writer.println(MODEL_HEADER);
			for(final Model model:models){
				for(int end=0; end<2; end++){
					final IntList list=model.errors[end];
					buffer.clear().append(members[model.member]).tab().append(model.name).tab().append(end==0 ? "5prime" : "3prime").tab().append(list.size);
					if(list.size==0){buffer.append("\tNA\tNA\tNA\tNA\t0\tNA\n");}
					else{
						final int q25=quantile(list, 25);
						final int q75=quantile(list, 75);
						int within=0;
						for(int i=0; i<list.size; i++){if(Math.abs((long)list.array[i]-model.median[end])<=3){within++;}}
						buffer.tab().append(model.median[end]).tab().append(q25).tab().append(q75).tab().append((long)q75-q25).tab().append(within).tab();
						percentage(buffer, within, list.size).nl();
					}
					output.writer.print(buffer);
				}
			}
		}
	}

	/** Catalog-wide support includes models that never won a retained match. */
	private void writeSupport(){
		try(Output output=new Output(out+".support.tsv")){
			output.writer.println("cohort\tmodels\tn0\tn1_4\tn5_9\tn10_19\tn20_49\tn50_99\tn100_plus\tminimum_n\tmedian_n\tmaximum_n");
			for(int member=0; member<=members.length; member++){
				final IntList counts=new IntList();
				final int[] bins=new int[7];
				for(final Model model:models){
					if(member!=members.length && model.member!=member){continue;}
					final int n=model.errors[0].size;
					counts.add(n);
					bins[n==0 ? 0 : n<5 ? 1 : n<10 ? 2 : n<20 ? 3 : n<50 ? 4 : n<100 ? 5 : 6]++;
				}
				counts.sort();
				buffer.clear().append(label(member)).tab().append(counts.size);
				for(final int value:bins){buffer.tab().append(value);}
				if(counts.size==0){buffer.append("\tNA\tNA\tNA\n");}
				else{buffer.tab().append(counts.array[0]).tab().append(quantile(counts, 50)).tab().append(counts.array[counts.size-1]).nl();}
				output.writer.print(buffer);
			}
		}
	}

	/** Describes median-size/spread strata without selecting a new shifting rule.
	 * Single-end bins use that end; paired bins use the maximum absolute applied
	 * median and maximum IQR across the two ends. Unsupported/invalid spans apply
	 * zero shift. Training IQR is still reported separately, including low support.
	 */
	private void writeBreakdown(){
		final Stats[][][][][] bins=new Stats[panels.length][3][2][5][6];
		for(int panel=0; panel<panels.length; panel++){
			long total=0;
			long selected=0;
			final ObjectSet<Target> seenReferences=new ObjectSet<Target>(Target.class);
			final ObjectSet<Target> seenCalls=new ObjectSet<Target>(Target.class);
			try(Input input=new Input(evaluation[panel])){
				RrnaResourceIO.header(input.reader, HEADER, evaluation[panel]);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					readRow(line, seenReferences, seenCalls);
					total++;
					final Model model=models.get(modelIndex);
					if(model.member!=breakdownIndex){continue;}
					selected++;
					final int support=model.errors[0].size>=thresholds[0] ? 1 : 0;
					final long m5=support==1 ? model.median[0] : 0;
					final long m3=support==1 ? model.median[1] : 0;
					final long shiftedStart=(long)callStart+(minus ? m3 : -m5);
					final long shiftedStop=(long)callStop+(minus ? m5 : -m3);
					final boolean held=shiftedStart<1 || shiftedStop<shiftedStart || shiftedStop>Integer.MAX_VALUE;
					final long a5=held ? error5 : error5-m5;
					final long a3=held ? error3 : error3-m3;
					final long iqr5=model.iqr(0);
					final long iqr3=model.iqr(1);
					for(int end=0; end<3; end++){
						final long magnitude=held ? 0 : end==0 ? Math.abs(m5) : end==1 ? Math.abs(m3) : Math.max(Math.abs(m5), Math.abs(m3));
						final long spread=end==0 ? iqr5 : end==1 ? iqr3 : Math.max(iqr5, iqr3);
						final int medianBin=bin(magnitude);
						final int spreadBin=spread<0 ? 5 : bin(spread);
						Stats stats=bins[panel][end][support][medianBin][spreadBin];
						if(stats==null){bins[panel][end][support][medianBin][spreadBin]=stats=new Stats();}
						if(end==0){stats.accept(error5, 0, a5, 0, support==1, held);}
						else if(end==1){stats.accept(error3, 0, a3, 0, support==1, held);}
						else{stats.accept(error5, error3, a5, a3, support==1, held);}
					}
				}
			}
			require(total==expected[panel+2] && selected==breakdownCounts[panel],
				"Breakdown must retain the accepted complete and member populations: "+panels[panel]+" "+total+"/"+selected);
		}
		try(Output output=new Output(out+".breakdown.tsv")){
			output.writer.println("support_group\tabs_applied_median_bin_bp\ttraining_iqr_bin_bp\t"+RESULT_HEADER);
			for(int panel=0; panel<panels.length; panel++){
				for(int end=0; end<3; end++){
					long conserved=0;
					for(int support=0; support<2; support++){
						for(int median=0; median<5; median++){
							for(int spread=0; spread<6; spread++){
								final Stats stats=bins[panel][end][support][median][spread];
								if(stats==null){continue;}
								conserved+=stats.n;
								buffer.clear().append(support==1 ? "SUPPORTED" : "UNSUPPORTED").tab().append(BINS[median]).tab().append(BINS[spread]).tab();
								stats.row(buffer, panels[panel], breakdownMember, thresholds[0], end);
								output.writer.print(buffer);
							}
						}
					}
					require(conserved==breakdownCounts[panel], "Median/spread bins must partition every fixed matched member row exactly once per end");
				}
			}
		}
	}

	/** Replays only training calls after learning medians; outer assemblies never
	 * enter residual histograms or weights. A second small witness pass retains
	 * the real coordinate geometry, which sorted per-model errors alone cannot do.
	 */
	private void writeCenteredHistograms(final Sources sources){
		final IntLongHashMap[][] histograms=new IntLongHashMap[members.length][2];
		final long[] retained=new long[members.length];
		final long[] supported=new long[members.length];
		final long[] held=new long[members.length];
		final long[] coordinateHeld=new long[members.length];
		final long[] belowMinimum=new long[members.length];
		final long[] rawBelowMinimum=new long[members.length];
		for(int member=0; member<members.length; member++){
			for(int end=0; end<2; end++){histograms[member][end]=new IntLongHashMap();}
		}
		final ObjectSet<Target> references=new ObjectSet<Target>(Target.class);
		final ObjectSet<Target> calls=new ObjectSet<Target>(Target.class);
		long selected=0;
		long included=0;
		for(final String path:training){
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, HEADER, path);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					readRow(line, references, calls);
					selected++;
					if(sources.excludes(row)){continue;}
					included++;
					final Model model=models.get(modelIndex);
					final boolean hasSupport=model.errors[0].size>=thresholds[0];
					final long m5=hasSupport ? model.median[0] : 0;
					final long m3=hasSupport ? model.median[1] : 0;
					final long start=(long)callStart+(minus ? m3 : -m5);
					final long stop=(long)callStop+(minus ? m5 : -m3);
					final boolean invalidCoordinates=start<1 || stop<start || stop>Integer.MAX_VALUE;
					final long minimum=minimumLengths==null ? 1 : minimumLengths[model.member];
					final boolean shortSpan=!invalidCoordinates && stop-start+1<minimum;
					final boolean invalid=invalidCoordinates || shortSpan;
					final long after5=invalid ? error5 : error5-m5;
					final long after3=invalid ? error3 : error3-m3;
					require(after5>=-Integer.MAX_VALUE && after5<=Integer.MAX_VALUE
						&& after3>=-Integer.MAX_VALUE && after3<=Integer.MAX_VALUE,
						"Valid one-based shifted coordinates bound signed errors for the existing integer histogram format");
					histograms[model.member][0].add((int)after5);
					histograms[model.member][1].add((int)after3);
					retained[model.member]++;
					if(hasSupport){supported[model.member]++;}
					if(invalid){held[model.member]++;}
					if(invalidCoordinates){coordinateHeld[model.member]++;}
					if(shortSpan){belowMinimum[model.member]++;}
					if((long)callStop-callStart+1<minimum){rawBelowMinimum[model.member]++;}
				}
			}
		}
		require(selected==expected[0] && included==expected[1], "The second witness pass must retain the same training-only population");
		for(int member=0; member<members.length; member++){
			require(held[member]==coordinateHeld[member]+belowMinimum[member],
				"Coordinate and minimum-length holds must partition all retained raw-error substitutions");
			for(int end=0; end<2; end++){
				require(histograms[member][end].sum()==retained[member], "Every retained member match must contribute exactly one residual per end");
			}
		}
		for(int member=0; member<members.length; member++){
			try(Output output=new Output(out+"."+members[member]+".histogram.tsv")){
				output.writer.println("family\tend\tmeasure\toffset_bp\tloci");
				for(int end=0; end<2; end++){
					final IntLongHashMap histogram=histograms[member][end];
					final int[] offsets=histogram.toArray();
					java.util.Arrays.sort(offsets);
					for(final int offset:offsets){
						output.writer.print(buffer.clear().append(members[member]).tab().append(end==0 ? "5prime" : "3prime")
							.append("\tsigned\t").append(offset).tab().append(histogram.get(offset)).nl());
					}
				}
			}
		}
		try(Output offsets=new Output(out+".offsets.tsv"); Output accounting=new Output(out+".training.tsv")){
			offsets.writer.println("member\tmodel\ttraining_n\tminimum_n\tapplied_median5_bp\tapplied_median3_bp");
			for(final Model model:models){
				final boolean active=model.errors[0].size>=thresholds[0];
				offsets.writer.print(buffer.clear().append(members[model.member]).tab().append(model.name)
					.tab().append(model.errors[0].size).tab().append(thresholds[0])
					.tab().append(active ? model.median[0] : 0).tab().append(active ? model.median[1] : 0).nl());
			}
			buffer.clear().append("member\ttraining_retained\tsupported\tgeometry_held\tminimum_n");
			if(minimumLengths!=null){buffer.append("\tminimum_length\tcoordinate_held\tbelow_minimum_held\traw_below_minimum");}
			accounting.writer.print(buffer.nl());
			for(int member=0; member<members.length; member++){
				buffer.clear().append(members[member]).tab().append(retained[member])
					.tab().append(supported[member]).tab().append(held[member]).tab().append(thresholds[0]);
				if(minimumLengths!=null){
					buffer.tab().append(minimumLengths[member]).tab().append(coordinateHeld[member])
						.tab().append(belowMinimum[member]).tab().append(rawBelowMinimum[member]);
				}
				accounting.writer.print(buffer.nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Disjoint change bins are meaningful only where both medians meet support. */
	private static final class Stability{
		void add(final boolean oldSupport, final boolean newSupport, final long delta){
			assert(!oldSupport || newSupport) : "Nested pooled training populations cannot lose supported models";
			models++;
			if(oldSupport){before++;}
			if(newSupport){after++;}
			if(!oldSupport && newSupport){gained++;}
			if(oldSupport && newSupport){
				both++;
				final long change=Math.abs(delta);
				if(change==0){same++;}
				else if(change<=2){small++;}
				else if(change<=5){moderate++;}
				else{large++;}
			}
			assert(same+small+moderate+large==both && before==both && after==both+gained)
				: "Stability bins and newly supported models must partition the nested support populations";
		}
		long models, before, after, both, gained, same, small, moderate, large;
	}

	private static final class Model{
		Model(final String name_, final int member_){
			name=name_;
			member=member_;
		}
		void finish(){
			assert(errors[0].size==errors[1].size) : "Each retained matched call contributes exactly one error per end";
			for(int end=0; end<2; end++){
				errors[end].sort();
				median[end]=errors[end].size==0 ? 0 : quantile(errors[end], 50);
			}
		}
		long iqr(final int end){
			return errors[end].size==0 ? -1 : (long)quantile(errors[end], 75)-quantile(errors[end], 25);
		}
		final String name;
		final int member;
		final IntList[] errors={new IntList(), new IntList()};
		final int[] median=new int[2];
	}

	/** Aggregates one end, or both ends using pooled RMS and summed squared distance. */
	private static final class Stats{
		void accept(final long b5, final long b3, final long a5, final long a3, final boolean hasSupport, final boolean held){
			assert(b5>=-Integer.MAX_VALUE && b5<=Integer.MAX_VALUE && b3>=-Integer.MAX_VALUE && b3<=Integer.MAX_VALUE
				&& a5>=-Integer.MAX_VALUE && a5<=Integer.MAX_VALUE && a3>=-Integer.MAX_VALUE && a3<=Integer.MAX_VALUE)
				: "Validated one-based int coordinates, including shifted-span holds, bound each error so two squared errors fit in long";
			final long before=b5*b5+b3*b3;
			final long after=a5*a5+a3*a3;
			final long b=Math.max(Math.abs(b5), Math.abs(b3));
			final long a=Math.max(Math.abs(a5), Math.abs(a3));
			n++;
			if(hasSupport){supported++;}
			if(b5!=a5 || b3!=a3){shifted++;}
			if(held){geometryHeld++;}
			if(before>after){better++;}
			else if(before<after){worse++;}
			else{tied++;}
			if(b==0){exactBefore++;}
			if(a==0){exactAfter++;}
			if(b<=3){within3Before++;}
			if(a<=3){within3After++;}
			if(b<=10){within10Before++;}
			if(a<=10){within10After++;}
			squaredBefore+=before;
			squaredAfter+=after;
			assert(better+worse+tied==n && exactAfter<=within3After && within3After<=within10After && within10After<=n)
				: "Disjoint movement outcomes and nested tolerances must conserve the fixed matched cohort";
		}
		void row(final ByteBuilder text, final String panel, final String cohort, final int threshold, final int end){
			text.append(panel).tab().append(cohort).tab().append(threshold).tab().append(end==2 ? "both" : end==0 ? "5prime" : "3prime")
				.tab().append(n).tab().append(supported).tab().append(shifted).tab().append(geometryHeld);
			for(final long value:new long[]{exactBefore, exactAfter, within3Before, within3After, within10Before, within10After}){
				text.tab().append(value).tab();
				percentage(text, value, n);
			}
			if(n==0){text.append("\tNA\tNA\tNA");}
			else{
				final double before=Math.sqrt(squaredBefore/(n*(end==2 ? 2.0 : 1.0)));
				final double after=Math.sqrt(squaredAfter/(n*(end==2 ? 2.0 : 1.0)));
				text.tab().append(before, 8).tab().append(after, 8).tab().append(after-before, 8);
			}
			text.tab().append(better).tab().append(worse).tab().append(tied).tab();
			if(better+worse==0){text.append("NA");}
			else{text.append(better/(double)(better+worse), 8);}
			text.tab().append(QUALIFICATION).nl();
		}
		long n, supported, shifted, geometryHeld, exactBefore, exactAfter, within3Before, within3After, within10Before, within10After;
		long better, worse, tied;
		double squaredBefore, squaredAfter;
	}

	/** Exact assembly exclusions use input bindings, never inferred TaxIDs or names. */
	private static final class Sources{
		Sources(final String bindings, final String audit){
			final ArrayList<String> assemblies=new ArrayList<String>();
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(bindings)){
				RrnaResourceIO.header(input.reader, "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff", bindings);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					require(p.terms()==5 && p.parseInt(0)==assemblies.size() && p.length(1)>0 && !p.termEquals(".", 1), "Source ordinals/assembly identities must be explicit and dense");
					probe.set(p, 2);
					require(p.length(2)>0 && !byInput.contains(probe), "Input IDs must bind uniquely to assemblies");
					byInput.put(probe.copy(), assemblies.size());
					assemblies.add(p.parseString(1));
				}
			}
			require(!assemblies.isEmpty(), "Offset learning requires a complete nonempty source map even without held-out assemblies");
			excluded=new boolean[assemblies.size()];
			if(audit==null){return;}// Only the explicit all-data refit constructor permits an absent audit.
			final ObjectSet<String> omitted=new ObjectSet<String>(String.class);
			final ObjectIntMap<String> panelRows=new ObjectIntMap<String>(String.class);
			try(Input input=new Input(audit)){
				RrnaResourceIO.header(input.reader, NcrnaEndpointFoldInputs.AUDIT_LABEL, audit);
				RrnaResourceIO.header(input.reader, NcrnaEndpointFoldInputs.AUDIT_HEADER, audit);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					require(p.terms()==5 && p.termEquals("EXCLUDED_ASSEMBLY", 4) && p.termStartsWith("/", 0) && p.termStartsWith("/", 2)
						&& p.length(3)>0 && !p.termEquals(".", 3) && !p.termEquals("NA", 3),
						"Every outer-panel row must bind an exact excluded assembly and absolute paths");
					final String panel=p.parseString(0);
					final int previous=panelRows.get(panel);
					require(p.parseInt(1)==(previous<0 ? 1 : previous+1), "Outer-panel rows must be dense and ordered");
					panelRows.put(panel, previous<0 ? 1 : previous+1);
					omitted.add(p.parseString(3));
				}
			}
			require(!assemblies.isEmpty() && panelRows.size()==2 && !omitted.isEmpty(), "Complete source map and both outer-panel exclusions are required");
			for(int i=0; i<excluded.length; i++){excluded[i]=omitted.contains(assemblies.get(i));}
		}
		boolean excludes(final LineParser1 row){
			probe.set(row, 0);
			final int index=byInput.get(probe);
			require(index>=0, "Training witness source is absent from its accepted assembly binding");
			return excluded[index];
		}
		final Target probe=new Target();
		final ObjectIntMap<Target> byInput=new ObjectIntMap<Target>(Target.class);
		final boolean[] excluded;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Nearest-rank quantiles of signed errors; IntList.percentile is sum-weighted. */
	private static int quantile(final IntList sorted, final int percent){
		assert(sorted.size>0 && percent>0 && percent<=100) : "Nearest-rank lookup requires a nonempty sorted population";
		return sorted.array[(int)(((long)sorted.size*percent+99)/100)-1];
	}
	private static int bin(final long size){
		assert(size>=0) : "Magnitude and sorted interquartile spread cannot be negative";
		return size==0 ? 0 : size<=2 ? 1 : size<=5 ? 2 : size<=10 ? 3 : 4;
	}
	private static ByteBuilder field(final ByteBuilder text, final LineParser1 row, final int index){
		row.setBounds(index);
		return text.append(row.line(), row.a(), row.b()-row.a());
	}
	private static ByteBuilder percentage(final ByteBuilder text, final long value, final long total){
		return total==0 ? text.append("NA") : text.append(100.0*value/total, 8);
	}
	private static long[] counts(final String text){
		final String[] fields=text.split(",", -1);
		final long[] values=new long[fields.length];
		for(int i=0; i<values.length; i++){
			values[i]=Long.parseLong(fields[i]);
			require(values[i]>=0, "Accepted counts cannot be negative");
		}
		return values;
	}
	private String label(final int member){return member==members.length ? family : members[member];}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String[] training, refs, members, evaluation, panels, comparison;
	private String family, bindings, exclude, out;
	private String breakdownMember;
	private long[] breakdownCounts;
	private long[] minimumLengths;
	private int breakdownIndex=-1;
	private boolean centeredHistograms, refit;
	private long[] expected;
	private int expectedModels;
	private int[] thresholds={5,10,20};
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ObjectIntMap<Target> byModel=new ObjectIntMap<Target>(Target.class);
	private final Target probe=new Target();
	private final LineParser1 row=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private int modelIndex, callStart, callStop, error5, error3;
	private boolean minus;
	private static final String HEADER="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
	private static final String MODEL_HEADER="member\tmodel\tend\tn\tmedian_error_bp\tq25_bp\tq75_bp\tiqr_bp\twithin3_of_median\twithin3_of_median_pct";
	private static final String QUALIFICATION="SAVED_MATCHES_NO_REMATCH_NO_CONTIG_LENGTHS_NO_PATH_SELECTION";
	private static final String[] BINS={"0", "1_2", "3_5", "6_10", "ABOVE_10", "NO_SUPPORT"};
	private static final String RESULT_HEADER="panel\tcohort\tminimum_model_n\tend\tmatched\tsupported\tshifted\tgeometry_held"
		+"\texact_before\texact_before_pct\texact_after\texact_after_pct\twithin3_before\twithin3_before_pct\twithin3_after\twithin3_after_pct"
		+"\twithin10_before\twithin10_before_pct\twithin10_after\twithin10_after_pct\trms_before_bp\trms_after_bp\trms_change_bp"
		+"\tbetter\tworse\ttied\tbetter_fraction\tqualification";
}
