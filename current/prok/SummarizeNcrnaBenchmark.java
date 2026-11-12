package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.IntHashMap;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/** Summarizes a declared ordinary-caller benchmark from saved native reports.
 * Derived from template.A_SampleByteFile. One thread owns the readers, interned
 * byte keys and primitive observations. Only distinct retained keys are copied.
 * No sequence reads, calling, alignment, regrading or truth admission occurs.
 *
 * Every declared policy/worker/repeat cell must occur once. Native summary counts
 * are reconciled with raw references, queries and input manifests before output.
 * Detection identity is the exact encoded input/family/reference triple; equal
 * counts do not imply equal sets. Compared cells must share the same declared
 * input bindings, reference parents, masks and reciprocal-overlap rule.
 * Explicit timingonly=t instead reads only timings and optional profile data,
 * rejects grade arguments and makes no accuracy or detection-stability claim.
 * Optional literal taxonomy metadata groups the same checked detections by its
 * declared domain/phylum/category/status, preserving zero-query inputs and unknown
 * labels. Taxonomy does not alter scoring or infer rarity or biological identity.
 * @author Brian Bushnell, Raiden, Ganyu
 */
public final class SummarizeNcrnaBenchmark{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final SummarizeNcrnaBenchmark tool=new SummarizeNcrnaBenchmark(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	public SummarizeNcrnaBenchmark(final String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("gradeinputs")){inputManifest=b;}
			else if(a.equals("timingonly")){timingOnly=Parse.parseBoolean(b);}
			else if(a.equals("grades")){summaries=b;}
			else if(a.equals("raw")){raw=b;}
			else if(a.equals("counts")){originalCounts=b;}
			else if(a.equals("taxonomy")){taxonomy=b;}
			else if(a.equals("subtyperaw")){subtypeRaw=b;}
			else if(a.equals("subtypegrades")){subtypeGrades=b;}
			else if(a.equals("profile")){profile=b;}
			else if(a.equals("profiletime")){profileTime=b;}
			else if(a.equals("baseline")){baseline=b;}
			else if(a.equals("policies")){policyNames=b==null ? null : b.split(",", -1);}
			else if(a.equals("workers")){workerCounts=Parse.parseIntArray(b, ",");}
			else if(a.equals("repeats")){repeats=Integer.parseInt(b);}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		timings=common.in1;
		prefix=common.out1;
		require(timings!=null && prefix!=null, "Supply in=timings out=fresh_prefix");
		if(timingOnly){
			require(subtypeRaw==null && subtypeGrades==null, "timingonly=t does not accept subtype reports");
			require(inputManifest==null && summaries==null && raw==null && originalCounts==null && taxonomy==null,
				"timingonly=t does not accept grade or count inputs or taxonomy; no accuracy or stability analysis is performed");
		}else{
			require(inputManifest!=null && summaries!=null && raw!=null, "Full analysis requires gradeinputs= grades= raw=");
		}
		require(policyNames!=null && policyNames.length>0 && repeats>0 && workerCounts!=null && workerCounts.length>0,
			"Declare policies=names workers=positive_counts repeats=positive_count");
		require(!common.overwrite && !common.append, "Benchmark reports require fresh outputs");
		require((profile==null)==(profileTime==null), "profile= and profiletime= must be supplied together");
		require((subtypeRaw==null)==(subtypeGrades==null), "subtyperaw= and subtypegrades= must be supplied together");
		policyBytes=new byte[policyNames.length][];
		int found=-1;
		for(int i=0; i<policyNames.length; i++){
			require(!policyNames[i].isEmpty(), "Policy names must be nonempty");
			require(policyNames[i].indexOf('\t')<0 && policyNames[i].indexOf('\r')<0 && policyNames[i].indexOf('\n')<0,
				"Policy names must fit one literal TSV field");
			policyBytes[i]=policyNames[i].getBytes(StandardCharsets.UTF_8);
			if(policyNames[i].equals(baseline)){found=i;}
			for(int j=0; j<i; j++){require(!policyNames[i].equals(policyNames[j]), "Duplicate declared policy");}
		}
		require(found>=0, "baseline= must name one declared policy");
		baselineIndex=found;
		for(int i=0; i<workerCounts.length; i++){
			require(workerCounts[i]>0, "Worker counts must be positive");
			for(int j=0; j<i; j++){require(workerCounts[i]!=workerCounts[j], "Duplicate declared worker count");}
		}
		cells=new Cell[Math.multiplyExact(Math.multiplyExact(policyNames.length, workerCounts.length), repeats)];
		for(int p=0; p<policyNames.length; p++){
			for(int t=0; t<workerCounts.length; t++){
				for(int r=1; r<=repeats; r++){cells[index(p, t, r)]=new Cell(p, t, r);}
			}
		}
		require(Tools.testInputFiles(false, true, timings, inputManifest, summaries, raw, originalCounts, taxonomy, profile, profileTime, subtypeRaw, subtypeGrades),
			"Every supplied benchmark input must be readable");
		if(timingOnly){RrnaResourceIO.fresh(prefix, ".cells.tsv", ".timing.tsv", ".profile.tsv");}
		else{RrnaResourceIO.fresh(prefix, ".cells.tsv", ".timing.tsv", ".grades.tsv", ".comparisons.tsv", ".changes.tsv", ".profile.tsv");}
		if(taxonomy!=null){RrnaResourceIO.fresh(prefix, ".lineages.tsv", ".input_changes.tsv");}
		if(subtypeRaw!=null){
			RrnaResourceIO.fresh(prefix, ".subtype_comparisons.tsv", ".subtype_changes.tsv");
			subtypeComparison=new NcrnaSubtypeComparison(cells.length);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Native readers        ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		readTimings();
		validateTiming();
		if(!timingOnly){
			readInputs();
			if(taxonomy!=null){readTaxonomy();}
			readSummaries();
			if(originalCounts!=null){readOriginalCounts();}
			readRaw();
			validate();
			if(subtypeComparison!=null){readSubtypeReports();}
		}
		writeTiming();
		if(!timingOnly){
			writeGrades();
			writeComparisons();
		}
		if(profile!=null){writeProfile();}
		timer.stop();
		if(timingOnly){
			outstream.println("Summarized timing for "+cells.length+" declared benchmark cells. Time: "+timer);
			outstream.println("Accuracy and detection stability were not evaluated; no grade inputs were read.");
			outstream.println("NCRNA_BENCHMARK_TIMING_ONLY_COMPLETE");
			return;
		}
		outstream.println("Summarized "+cells.length+" declared benchmark cells, "+inputs.size()+" inputs and "+loci.size()+" reference keys. Time: "+timer);
		outstream.println("Saved-report analysis only; no calling, regrading or biological reference decisions.");
		outstream.println("NCRNA_BENCHMARK_SUMMARY_COMPLETE: differing repeat comparisons="+unstableRepeats+
			"; differing worker comparisons="+unstableWorkers);
	}

	private void readTimings(){
		try(Input in=new Input(timings)){
			header(in, timings, TIMING_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 6);
				final Cell cell=cell();
				check(!cell.timed, "Each declared timing cell must occur exactly once");
				cell.user=number(3);
				cell.system=number(4);
				cell.wall=number(5);
				check(cell.user>=0 && cell.system>=0 && cell.wall>0 && finite(cell.user+cell.system) && cell.user+cell.system>0,
					"Timing requires nonnegative user/system and positive finite CPU/wall seconds");
				cell.timed=true;
			}
		}
		for(final Cell cell:cells){
			if(!cell.timed){check(false, "A declared timing cell is missing: "+cellLabel(cell));}
		}
	}

	private void readInputs(){
		try(Input in=new Input(inputManifest)){
			header(in, inputManifest, CONTEXT+GradeNcrnaReference.INPUT_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 13);
				final Cell cell=cell();
				final int id=inputs.intern(parser, 3);
				check(cell.inputRows.get(id)<0, "An input must appear once per benchmark cell");
				if(id==fastaPaths.size()){
					fastaPaths.add(copy(4));
					referencePaths.add(copy(6));
				}else{
					check(matches(fastaPaths.get(id), 4) && matches(referencePaths.get(id), 6),
						"Repeated inputs must retain the exact FASTA and reference-GFF bindings");
				}
				final long features=count(7), included=count(11), excluded=count(12);
				check(features==Math.addExact(included, excluded) && count(8)==0 && count(9)==0 && count(10)==0,
					"Bound query manifests must conserve their generic-RNA projection");
				cell.inputRows.put(id, cell.inputIds.size);
				cell.inputIds.add(id);
				cell.expectedQueries.add(features);
			}
		}
		for(final Cell cell:cells){
			check(cell.inputIds.size>0, "Every cell needs a complete declared input population, including empty queries");
			cell.actualQueries=new long[cell.inputIds.size];
			cell.sourceSeen=new boolean[cell.inputIds.size];
		}
	}

	/** Binds literal taxonomy fields to the native reports' escaped input/path keys.
	 * A larger metadata table is allowed; only selected inputs require one row.
	 * Category and unknown/status labels are supplied evidence, never inferred.
	 */
	private void readTaxonomy(){
		lineageByInput=new int[inputs.size()];
		Arrays.fill(lineageByInput, -1);
		corpusInputs=new long[inputs.size()];
		try(Input in=new Input(taxonomy)){
			header(in, taxonomy, TAXONOMY_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 7);
				nativeField(keyBuffer.clear(), 0);
				probe.set(keyBuffer.array, 0, keyBuffer.length());
				final int id=inputs.index.get(probe);
				if(id<0){continue;}
				check(lineageByInput[id]<0, "Each selected input needs exactly one taxonomy row");
				nativeField(keyBuffer.clear(), 1);
				probe.set(keyBuffer.array, 0, keyBuffer.length());
				check(fastaPaths.get(id).equals(probe), "Taxonomy must retain the selected input's exact FASTA path");
				corpusInputs[id]=count(4);
				keyBuffer.clear();
				for(final int field:TAXONOMY_LABEL_FIELDS){
					parser.setBounds(field);
					check(parser.b()>parser.a(), "Taxonomy labels must be explicit, including unknown labels");
					if(keyBuffer.length()>0){keyBuffer.tab();}
					nativeField(keyBuffer, field);
				}
				probe.set(keyBuffer.array, 0, keyBuffer.length());
				lineageByInput[id]=lineages.intern(probe);
			}
		}
		lineagePopulation=new int[lineages.size()];
		for(final int group:lineageByInput){
			check(group>=0, "Every selected input needs declared taxonomy, including zero-query inputs");
			lineagePopulation[group]++;
		}
	}

	/** Reuses the grade producer's byte encoder; no per-field String is created. */
	private ByteBuilder nativeField(final ByteBuilder out, final int field){
		parser.setBounds(field);
		return BuildNcrnaReference.field(out, parser.line(), parser.a(), parser.b());
	}

	private void readSummaries(){
		try(Input in=new Input(summaries)){
			header(in, summaries, CONTEXT+GradeNcrnaReference.SUMMARY_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 20);
				final Cell cell=cell();
				final int family=families.intern(parser, 3);
				while(cell.grades.size()<=family){cell.grades.add(null);}
				check(cell.grades.get(family)==null, "Each cell needs one summary per family");
				final Grade grade=new Grade();
				for(int i=0; i<GRADE_COUNTS; i++){grade.expected[i]=count(i+4);}
				check(grade.expected[TP]+grade.expected[FN]==grade.expected[SCORED], "Native summary TP+FN must equal scored references");
				check(grade.expected[BOTH]<=grade.expected[EXACT5] && grade.expected[BOTH]<=grade.expected[EXACT3]
					&& grade.expected[EXACT5]<=grade.expected[TP] && grade.expected[EXACT3]<=grade.expected[TP],
					"Exact endpoint counts must refer only to matched TPs");
				checkRatio(16, grade.expected[TP], Math.addExact(grade.expected[TP], grade.expected[FP]));
				checkRatio(17, grade.expected[TP], grade.expected[SCORED]);
				final double value=number(18);
				check(value>0 && value<=1 && (overlap<0 || Double.compare(overlap, value)==0),
					"All grades must use the same positive reciprocal-overlap rule");
				overlap=value;
				cell.grades.set(family, grade);
			}
		}
		for(final Cell cell:cells){
			check(cell.grades.size()==families.size(), "Every cell must retain every declared grade family, including exclusions");
			for(final Grade grade:cell.grades){check(grade!=null, "A cell is missing a family summary");}
		}
	}

	private void readOriginalCounts(){
		try(Input in=new Input(originalCounts)){
			header(in, originalCounts, CONTEXT+ORIGINAL_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 11);
				final Cell cell=cell();
				final int id=inputs.find(parser, 3), local=id<0 ? -1 : cell.inputRows.get(id);
				check(local>=0 && !cell.sourceSeen[local], "Original feature counts require each declared input once");
				long sum=0;
				for(int i=0; i<4; i++){
					final long value=count(5+i);
					cell.sourceTotals[i]=Math.addExact(cell.sourceTotals[i], value);
					sum=Math.addExact(sum, value);
				}
				check(sum==count(4) && count(8)==count(9) && count(9)==cell.expectedQueries.get(local),
					"Original feature types and projected queries must conserve the bound source population");
				cell.sourceSeen[local]=true;
			}
		}
	}

	private void readRaw(){
		try(Input in=new Input(raw)){
			header(in, raw, CONTEXT+GradeNcrnaReference.RAW_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 21);
				final Cell cell=cell();
				if(subtypeComparison!=null){
					try{subtypeComparison.ordinary(index(cell.policy, cell.worker, cell.repeat), parser, 3);}
					catch(IllegalArgumentException failure){throw subtypeFailure(failure);}
				}
				final int input=inputs.find(parser, 3), family=families.find(parser, 5);
				final int local=input<0 ? -1 : cell.inputRows.get(input);
				check(local>=0 && family>=0, "Raw rows require an input and family declared by this cell");
				final Grade grade=cell.grades.get(family);
				if(parser.termEquals("REFERENCE", 4)){
					final int id=parent(input, family, 11);
					cell.references.add(id);
					final boolean masked=parser.termEquals("DONT_CARE", 10);
					final Locus locus=loci.get(id);
					check(locus.status<0 || locus.status==(masked ? 1 : 0), "A reference parent cannot change scored/masked status between cells");
					locus.status=masked ? 1 : 0;
					if(masked){grade.observed[MASKED_REF]++;}
					else{
						grade.observed[SCORED]++;
						if(parser.termEquals("TP", 10)){
							grade.observed[TP]++;
							cell.tps.add(id);
						}else{
							check(parser.termEquals("FN", 10), "A scored native reference outcome must be TP or FN");
							grade.observed[FN]++;
						}
					}
				}else{
					check(parser.termEquals("QUERY", 4), "Raw native rows must be REFERENCE or QUERY");
					cell.actualQueries[local]++;
					if(parser.termEquals("EXCLUDED", 10)){grade.observed[EXCLUDED]++;}
					else{
						if(parser.termEquals("true", 16)){grade.observed[DUPLICATES]++;}
						else{check(parser.termEquals("false", 16), "Evaluated queries require a native duplicate flag");}
						if(parser.termEquals("TP", 10)){
							check(!parser.termEquals(".", 13), "A query TP must name its scored reference parent");
							cell.queryTps.add(parent(input, family, 13));
							final boolean exact5=parser.parseInt(19)==0, exact3=parser.parseInt(20)==0;
							if(exact5){grade.observed[EXACT5]++;}
							if(exact3){grade.observed[EXACT3]++;}
							if(exact5 && exact3){grade.observed[BOTH]++;}
						}else if(parser.termEquals("FP", 10)){grade.observed[FP]++;}
						else{
							check(parser.termEquals("MASKED", 10), "Evaluated native query outcome must be TP, FP or MASKED");
							grade.observed[MASKED_CALL]++;
						}
					}
				}
			}
		}
	}

	/** Context framing matches the ordinary native tables; no new grading occurs. */
	private void readSubtypeReports(){
		assert(subtypeComparison!=null) : "Paired optional subtype paths create their checked saved-table consumer";
		try(Input in=new Input(subtypeGrades)){
			header(in, subtypeGrades, CONTEXT+NcrnaSubtypeGrade.SUMMARY_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 20);
				final Cell cell=cell();
				subtypeComparison.summary(index(cell.policy, cell.worker, cell.repeat), parser, 3, cell.inputIds.size, overlap);
			}
		}catch(IllegalArgumentException failure){throw subtypeFailure(failure);}
		try(Input in=new Input(subtypeRaw)){
			header(in, subtypeRaw, CONTEXT+NcrnaSubtypeGrade.RAW_HEADER);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				fields(row, 25);
				final Cell cell=cell();
				final int input=inputs.find(parser, 5);
				check(input>=0 && cell.inputRows.get(input)>=0, "Subtype row requires an input declared in this cell");
				subtypeComparison.raw(index(cell.policy, cell.worker, cell.repeat), parser, 3);
			}
		}catch(IllegalArgumentException failure){throw subtypeFailure(failure);}
		subtypeComparison.validate();
	}

	/** Keep the failing native path/line when the optional consumer rejects a row. */
	private IllegalArgumentException subtypeFailure(final IllegalArgumentException failure){
		assert(failure!=null) : "Only a caught subtype parsing failure is wrapped with its source location";
		return new IllegalArgumentException(failure.getMessage()+" at "+currentPath+":"+currentLine, failure);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Checked comparisons      ----------------*/
	/*--------------------------------------------------------------*/

	private void validate(){
		final Cell base=cells[index(baselineIndex, 0, 1)];
		for(final Cell cell:cells){
			check(cell.inputIds.size==base.inputIds.size, "Compared cells must contain the same declared inputs");
			for(int i=0; i<base.inputIds.size; i++){
				check(cell.inputRows.get(base.inputIds.get(i))>=0, "Compared cells must retain zero-query inputs too");
			}
			for(int i=0; i<cell.inputIds.size; i++){
				check(cell.actualQueries[i]==cell.expectedQueries.get(i), "Every bound generic-RNA query must occur once in the raw report");
				check(originalCounts==null || cell.sourceSeen[i], "Original feature counts are missing an input");
			}
			unique(cell.references);
			unique(cell.tps);
			unique(cell.queryTps);
			check(equal(cell.tps, cell.queryTps), "Native reference and query TPs must form the same one-to-one parent set");
			for(int family=0; family<cell.grades.size(); family++){
				final Grade grade=cell.grades.get(family);
				grade.observed[INPUTS]=cell.inputIds.size;
				if(!Arrays.equals(grade.expected, grade.observed)){
					check(false, "Native summary and raw counts disagree: "+cellLabel(cell)+" family="+
						new String(families.key(family).bytes, StandardCharsets.UTF_8)+" summary="+
						Arrays.toString(grade.expected)+" raw="+Arrays.toString(grade.observed));
				}
			}
		}
		for(final Cell cell:cells){
			check(equal(base.references, cell.references), "Compared cells must share the exact scored and masked reference-parent universe");
		}
	}

	/** The cluster timing floor applies whether or not saved grades are supplied. */
	private void validateTiming(){
		for(final Cell cell:cells){
			if(cell.policy==baselineIndex && workerCounts[cell.worker]==1){
				check(cell.user+cell.system>=10, "The shipped one-worker timing baseline must contain at least ten CPU seconds");
			}
		}
	}

	/** Merges sorted exact-parent IDs; first/last detection is counted per input and family. */
	private boolean compare(final String kind, final Cell left, final Cell right, final Output summary, final Output changes,
		final Output lineageOut, final Output inputChanges, final Output subtypeSummary, final Output subtypeChanges, final ComparisonScratch scratch){
		assert(left.tps.size==0 || left.tps.get(0)>=0) : "validate() establishes sorted nonnegative parent IDs for the comparison merge";
		final long[] lost=scratch.lost, gained=scratch.gained;
		Arrays.fill(lost, 0);
		Arrays.fill(gained, 0);
		Arrays.fill(scratch.lineageLost, 0);
		Arrays.fill(scratch.lineageGained, 0);
		int a=0, b=0;
		boolean different=false;
		while(a<left.tps.size || b<right.tps.size){
			final int x=a<left.tps.size ? left.tps.get(a) : Integer.MAX_VALUE;
			final int y=b<right.tps.size ? right.tps.get(b) : Integer.MAX_VALUE;
			if(x==y){a++; b++;}
			else{
				different=true;
				final boolean losing=x<y;
				final int id=losing ? left.tps.get(a++) : right.tps.get(b++);
				final Locus locus=loci.get(id);
				if(losing){lost[locus.family]++;}
				else{gained[locus.family]++;}
				if(taxonomy!=null){
					final int group=locus.family*lineages.size()+lineageByInput[locus.input];
					if(losing){scratch.lineageLost[group]++;}
					else{scratch.lineageGained[group]++;}
				}
				comparisonPrefix(buffer.clear(), kind, left, right).tab().append(losing ? "LOST" : "GAINED").tab()
					.append(parents.key(id).bytes).nl();
				changes.writer.print(buffer);
			}
		}
		final int[] before=scratch.before, after=scratch.after;
		for(int family=0; family<families.size(); family++){
			Arrays.fill(before, 0);
			Arrays.fill(after, 0);
			countParents(left, family, before);
			countParents(right, family, after);
			Arrays.fill(scratch.lineageCounts, 0);
			long first=0, last=0;
			for(int i=0; i<before.length; i++){
				final boolean newInput=before[i]==0 && after[i]>0, lostInput=before[i]>0 && after[i]==0;
				if(newInput){first++;}
				if(lostInput){last++;}
				if(taxonomy!=null){
					final int group=lineageByInput[i], offset=group*LINEAGE_COUNTS;
					final long[] counts=scratch.lineageCounts;
					counts[offset]+=before[i];
					counts[offset+1]+=after[i];
					if(before[i]>0){counts[offset+2]++;}
					if(after[i]>0){counts[offset+3]++;}
					if(newInput){counts[offset+4]++;}
					if(lostInput){counts[offset+5]++;}
					if(newInput || lostInput){
						comparisonPrefix(buffer.clear(), kind, left, right).tab().append(families.key(family).bytes).tab()
							.append(newInput ? "FIRST_DETECTION" : "LAST_LOSS").tab().append(inputs.key(i).bytes).tab()
							.append(lineages.key(group).bytes).tab().append(corpusInputs[i]).tab()
							.append(before[i]).tab().append(after[i]).nl();
						inputChanges.writer.print(buffer);
					}
				}
			}
			final long leftTp=left.grades.get(family).expected[TP], rightTp=right.grades.get(family).expected[TP];
			check(leftTp-lost[family]+gained[family]==rightTp, "Exact lost/gained parent counts must reconcile with native TPs");
			comparisonPrefix(buffer.clear(), kind, left, right).tab().append(families.key(family).bytes).tab()
				.append(leftTp).tab().append(rightTp).tab().append(lost[family]).tab().append(gained[family])
				.tab().append(first).tab().append(last).nl();
			summary.writer.print(buffer);
			if(taxonomy!=null){writeLineages(kind, left, right, family, lineageOut, scratch, first, last);}
		}
		if(subtypeComparison!=null){
			subtypeComparison.compare(comparisonPrefix(buffer.clear(), kind, left, right),
				index(left.policy, left.worker, left.repeat), index(right.policy, right.worker, right.repeat), subtypeSummary, subtypeChanges);
		}
		return different;
	}

	/** Each input belongs to one declared class; summed class counts must conserve
	 * the exact-parent comparison and its first/last input counts.
	 */
	private void writeLineages(final String kind, final Cell left, final Cell right, final int family,
		final Output out, final ComparisonScratch scratch, final long first, final long last){
		assert(lineageByInput.length==inputs.size()) : "readTaxonomy binds every input before class counts are emitted";
		long leftSum=0, rightSum=0, lostSum=0, gainedSum=0, firstSum=0, lastSum=0;
		for(int group=0; group<lineages.size(); group++){
			final int offset=group*LINEAGE_COUNTS, change=family*lineages.size()+group;
			final long[] counts=scratch.lineageCounts;
			final long lost=scratch.lineageLost[change], gained=scratch.lineageGained[change];
			check(counts[offset]-lost+gained==counts[offset+1], "Each taxonomy class must conserve its exact lost/gained references");
			leftSum+=counts[offset];
			rightSum+=counts[offset+1];
			lostSum+=lost;
			gainedSum+=gained;
			firstSum+=counts[offset+4];
			lastSum+=counts[offset+5];
			comparisonPrefix(buffer.clear(), kind, left, right).tab().append(families.key(family).bytes).tab()
				.append(lineages.key(group).bytes).tab().append(lineagePopulation[group]).tab()
				.append(counts[offset]).tab().append(counts[offset+1]).tab().append(lost).tab().append(gained);
			for(int i=2; i<LINEAGE_COUNTS; i++){buffer.tab().append(counts[offset+i]);}
			out.writer.print(buffer.nl());
		}
		check(leftSum==left.grades.get(family).expected[TP] && rightSum==right.grades.get(family).expected[TP]
			&& lostSum==scratch.lost[family] && gainedSum==scratch.gained[family] && firstSum==first && lastSum==last,
			"Declared taxonomy classes must conserve all family TP, change and first/last counts");
	}

	private void countParents(final Cell cell, final int family, final int[] counts){
		assert(counts.length==inputs.size()) : "Per-input detection counters must use the shared declared input dictionary";
		for(int i=0; i<cell.tps.size; i++){
			final Locus locus=loci.get(cell.tps.get(i));
			if(locus.family==family){counts[locus.input]++;}
		}
	}

	private void writeComparisons(){
		final ComparisonScratch scratch=new ComparisonScratch(inputs.size(), families.size(), lineages.size());
		try(Output summary=new Output(prefix+".comparisons.tsv"); Output changes=new Output(prefix+".changes.tsv");
			Output lineageOut=taxonomy==null ? null : new Output(prefix+".lineages.tsv");
			Output inputChanges=taxonomy==null ? null : new Output(prefix+".input_changes.tsv");
			Output subtypeSummary=subtypeComparison==null ? null : new Output(prefix+".subtype_comparisons.tsv");
			Output subtypeChanges=subtypeComparison==null ? null : new Output(prefix+".subtype_changes.tsv")){
			summary.writer.println(COMPARISON_HEADER+"\tfamily\tleft_TP\tright_TP\tlost_references\tgained_references\tfirst_detected_inputs\tlast_lost_inputs");
			changes.writer.println(COMPARISON_HEADER+"\tchange\tinput_id\tfamily\treference_id");
			if(subtypeComparison!=null){
				subtypeSummary.writer.println(COMPARISON_HEADER+"\t"+NcrnaSubtypeComparison.SUMMARY_TAIL);
				subtypeChanges.writer.println(COMPARISON_HEADER+"\t"+NcrnaSubtypeComparison.CHANGE_TAIL);
			}
			if(taxonomy!=null){
				lineageOut.writer.println(COMPARISON_HEADER+"\tfamily\tdomain\tphylum\tcategory\ttaxonomy_status\tselected_inputs"+
					"\tleft_TP\tright_TP\tlost_references\tgained_references\tleft_detected_inputs\tright_detected_inputs\tfirst_detected_inputs\tlast_lost_inputs");
				inputChanges.writer.println(COMPARISON_HEADER+"\tfamily\tchange\tinput_id\tdomain\tphylum\tcategory\ttaxonomy_status"+
					"\tdeclared_corpus_inputs\tleft_TP\tright_TP");
			}
			for(final Cell cell:cells){
				compare("POLICY_VS_BASELINE", cells[index(baselineIndex, cell.worker, cell.repeat)], cell, summary, changes, lineageOut, inputChanges, subtypeSummary, subtypeChanges, scratch);
				if(cell.repeat>1 && compare("REPEAT_VS_FIRST", cells[index(cell.policy, cell.worker, 1)], cell, summary, changes, lineageOut, inputChanges, subtypeSummary, subtypeChanges, scratch)){
					unstableRepeats++;
				}
				if(cell.worker>0 && compare("WORKERS_VS_FIRST", cells[index(cell.policy, 0, cell.repeat)], cell, summary, changes, lineageOut, inputChanges, subtypeSummary, subtypeChanges, scratch)){
					unstableWorkers++;
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Timing reports        ----------------*/
	/*--------------------------------------------------------------*/

	private void writeTiming(){
		final double[][][] stats=new double[policyNames.length*workerCounts.length][4][3];
		final double[] values=new double[repeats];
		for(int p=0; p<policyNames.length; p++){
			for(int t=0; t<workerCounts.length; t++){
				for(int metric=0; metric<4; metric++){
					for(int r=1; r<=repeats; r++){values[r-1]=seconds(cells[index(p, t, r)], metric);}
					Arrays.sort(values);
					final double[] result=stats[p*workerCounts.length+t][metric];
					result[0]=median(values);
					result[1]=values[0];
					result[2]=values[values.length-1];
				}
			}
		}
		try(Output rawCells=new Output(prefix+".cells.tsv"); Output reduced=new Output(prefix+".timing.tsv")){
			rawCells.writer.println("policy\tworkers\trepeat\tuser_seconds\tsystem_seconds\tcpu_seconds\twall_seconds\tCDS\ttRNA\trRNA\tgeneric_RNA");
			reduced.writer.print("policy\tworkers\trepeats\tuser_median\tuser_min\tuser_max\tsystem_median\tsystem_min\tsystem_max\tcpu_median\tcpu_min\tcpu_max\twall_median\twall_min\twall_max\tcpu_saved_percent\twall_saved_percent");
			if(timingOnly){reduced.writer.print("\tqualification");}
			reduced.writer.print("\n");
			for(final Cell cell:cells){
				cellPrefix(buffer.clear(), cell);
				for(int metric=0; metric<4; metric++){buffer.tab().append(seconds(cell, metric), 6);}
				for(final long value:cell.sourceTotals){
					buffer.tab();
					if(originalCounts==null){buffer.append("NA");}
					else{buffer.append(value);}
				}
				rawCells.writer.print(buffer.nl());
			}
			for(int p=0; p<policyNames.length; p++){
				for(int t=0; t<workerCounts.length; t++){
					buffer.clear().append(policyBytes[p]).tab().append(workerCounts[t]).tab().append(repeats);
					final double[][] result=stats[p*workerCounts.length+t];
					final double[][] base=stats[baselineIndex*workerCounts.length+t];
					for(final double[] metric:result){for(final double value:metric){buffer.tab().append(value, 6);}}
					buffer.tab().append(100*(1-result[2][0]/base[2][0]), 6)
						.tab().append(100*(1-result[3][0]/base[3][0]), 6);
					if(timingOnly){buffer.append("\tCaller outputs were not compared");}
					reduced.writer.print(buffer.nl());
				}
			}
		}
	}

	private void writeGrades(){
		try(Output out=new Output(prefix+".grades.tsv")){
			out.writer.println("policy\tworkers\trepeat\t"+GradeNcrnaReference.SUMMARY_HEADER);
			for(final Cell cell:cells){
				for(int f=0; f<families.size(); f++){
					final long[] counts=cell.grades.get(f).expected;
					cellPrefix(buffer.clear(), cell).tab().append(families.key(f).bytes);
					for(final long count:counts){buffer.tab().append(count);}
					buffer.tab();
					ratio(buffer, counts[TP], Math.addExact(counts[TP], counts[FP]));
					buffer.tab();
					ratio(buffer, counts[TP], counts[SCORED]);
					buffer.tab().append(overlap, 6).append("\tNATIVE_GRADE_RECONCILED_NOT_NEW_SCORING\n");
					out.writer.print(buffer);
				}
			}
		}
	}

	/** Worker stages are disjoint; the process denominator belongs to this profile run. */
	private void writeProfile(){
		final double processCpu;
		try(Input in=new Input(profileTime)){
			currentPath=profileTime;
			currentLine=0;
			final byte[] row=in.reader.nextLine();
			check(row!=null, "The profile needs one process user/system/wall timing row");
			fields(row, 3);
			final double user=number(0), system=number(1), wall=number(2);
			check(user>=0 && system>=0 && wall>0 && finite(user+system) && user+system>0,
				"Profile process timing needs positive finite CPU and elapsed seconds");
			processCpu=user+system;
			check(in.reader.nextLine()==null, "The profile denominator must describe exactly one invocation");
		}
		final Pool names=new Pool();
		final ArrayList<long[]> times=new ArrayList<long[]>();
		final IntList seen=new IntList();
		final long[] global=new long[5];
		int globalSeen=0;
		try(Input in=new Input(profile)){
			currentPath=profile;
			currentLine=0;
			byte[] row;
			do{
				row=in.reader.nextLine();
				currentLine++;
			}while(row!=null && row.length>0 && row[0]=='#');
			check(row!=null && Arrays.equals(row, PROFILE_HEADER.getBytes(StandardCharsets.US_ASCII)), "Unexpected native CPU-profile header");
			while((row=in.reader.nextLine())!=null){
				fields(row, 10);
				final long cpu=count(3);
				count(2);
				if(parser.termEquals("GLOBAL", 0)){
					int stage=-1;
					for(int i=0; i<GLOBAL_STAGES.length; i++){if(parser.termEquals(GLOBAL_STAGES[i], 1)){stage=i;}}
					check(stage>=0 && (globalSeen&(1<<stage))==0, "Each native global CPU stage must appear once");
					global[stage]=cpu;
					globalSeen|=1<<stage;
				}else{
					final int name=names.intern(parser, 0), stage=parser.termEquals("INIT", 1) ? 0 : parser.termEquals("SCAVENGE", 1) ? 1 : -1;
					if(name==times.size()){
						times.add(new long[2]);
						seen.add(0);
					}
					check(stage>=0 && (seen.get(name)&(1<<stage))==0, "Each named RNA family needs one INIT and SCAVENGE CPU span");
					times.get(name)[stage]=cpu;
					seen.set(name, seen.get(name)|(1<<stage));
				}
			}
		}
		check(globalSeen==31, "The native profile must include every global stage, including zero-work stages");
		try(Output out=new Output(prefix+".profile.tsv")){
			out.writer.println("component\tkind\tcpu_ns\tcpu_seconds\tprofile_process_cpu_seconds\tprocess_cpu_percent\tqualification");
			long named=0;
			for(int i=0; i<names.size(); i++){
				check(seen.get(i)==3, "A profiled family is missing INIT or SCAVENGE");
				final long sum=Math.addExact(times.get(i)[0], times.get(i)[1]);
				named=Math.addExact(named, sum);
				profileRow(out, names.key(i).bytes, "FAMILY_INIT_AND_SCAVENGE", sum, processCpu);
			}
			profileRow(out, ascii("NAMED_RNA_SUBTOTAL"), "SUBTOTAL_NOT_ADDITIONAL_WORK", named, processCpu);
			long covered=named;
			for(int i=0; i<global.length; i++){
				final boolean path=i==3;
				profileRow(out, ascii(GLOBAL_STAGES[i]), path ? "SHARED_CDS_RNA_EXCLUDED" : "SHARED_OR_EXTERNAL_RNA", global[i], processCpu);
				if(!path){covered=Math.addExact(covered, global[i]);}
			}
			profileRow(out, ascii("COVERED_RNA_TOTAL"), "COVERED_NOT_EXHAUSTIVE_ATTRIBUTION", covered, processCpu);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Parsing helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private void header(final Input in, final String path, final String expected){
		currentPath=path;
		currentLine=1;
		RrnaResourceIO.header(in.reader, expected, path);
	}

	private void fields(final byte[] row, final int expected){
		currentLine++;
		parser.set(row);
		check(parser.terms()==expected, "Unexpected native table field count");
	}

	private long count(final int field){
		final long value=parser.parseLong(field);
		check(value>=0, "Native observation counts and CPU nanoseconds must be nonnegative");
		return value;
	}

	private double number(final int field){
		final double value=parser.parseDouble(field);
		check(finite(value), "Reported numeric observations must be finite");
		return value;
	}

	/** The native GradeNcrnaReference writer displays ratios to six decimal places. */
	private void checkRatio(final int field, final long numerator, final long denominator){
		if(denominator==0){check(parser.termEquals("NA", field), "A zero native ratio denominator must be reported as NA");}
		else{
			check(Math.abs(number(field)-numerator/(double)denominator)<=0.000000500001,
				"Native precision/recall must agree with counts within six-place rounding");
		}
	}

	private Cell cell(){
		int policy=-1, worker=-1;
		for(int i=0; i<policyBytes.length; i++){
			probe.set(parser, 0);
			if(bytesEqual(probe, policyBytes[i])){
				policy=i;
				break;
			}
		}
		final int workers=parser.parseInt(1), repeat=parser.parseInt(2);
		for(int i=0; i<workerCounts.length; i++){if(workerCounts[i]==workers){worker=i;}}
		check(policy>=0 && worker>=0 && repeat>=1 && repeat<=repeats,
			"A row names a policy, worker count or repeat outside the declared benchmark roster");
		return cells[index(policy, worker, repeat)];
	}

	private int index(final int policy, final int worker, final int repeat){
		assert(policy>=0 && policy<policyNames.length && worker>=0 && worker<workerCounts.length && repeat>=1 && repeat<=repeats)
			: "The declared policy/worker/repeat product defines the cell array indices";
		return (policy*workerCounts.length+worker)*repeats+repeat-1;
	}

	/** Exact byte identity, not a digest: equal parent IDs in different inputs remain distinct. */
	private int parent(final int input, final int family, final int parentField){
		parser.setBounds(parentField);
		check(parser.b()>parser.a() && !parser.termEquals(".", parentField), "A native reference or matched query needs a stable parent ID");
		keyBuffer.clear();
		appendField(keyBuffer, 3).tab();
		appendField(keyBuffer, 5).tab();
		appendField(keyBuffer, parentField);
		probe.set(keyBuffer.array, 0, keyBuffer.length());
		final int id=parents.intern(probe);
		if(id==loci.size()){loci.add(new Locus(input, family));}
		assert(id<loci.size()) : "The parent dictionary and retained locus records grow together";
		return id;
	}

	private ByteBuilder appendField(final ByteBuilder out, final int field){
		parser.setBounds(field);
		return out.append(parser.line(), parser.a(), parser.b()-parser.a());
	}

	private Target copy(final int field){
		probe.set(parser, field);
		return probe.copy();
	}

	private boolean matches(final Target value, final int field){
		probe.set(parser, field);
		return value.equals(probe);
	}

	private void unique(final IntList values){
		values.sort();
		for(int i=1; i<values.size; i++){
			check(values.get(i-1)<values.get(i), "One cell repeats a native reference or consumes a scored parent more than once");
		}
	}

	private static boolean equal(final IntList left, final IntList right){
		if(left.size!=right.size){return false;}
		for(int i=0; i<left.size; i++){if(left.get(i)!=right.get(i)){return false;}}
		return true;
	}

	private static boolean bytesEqual(final Target key, final byte[] bytes){
		if(key.to-key.from!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(key.bytes[key.from+i]!=bytes[i]){return false;}}
		return true;
	}

	private ByteBuilder cellPrefix(final ByteBuilder out, final Cell cell){
		return out.append(policyBytes[cell.policy]).tab().append(workerCounts[cell.worker]).tab().append(cell.repeat);
	}

	/** Allocated only for a failure diagnostic, never for the record-processing path. */
	private String cellLabel(final Cell cell){
		return policyNames[cell.policy]+" workers="+workerCounts[cell.worker]+" repeat="+cell.repeat;
	}

	private ByteBuilder comparisonPrefix(final ByteBuilder out, final String kind, final Cell left, final Cell right){
		out.append(kind).tab();
		cellPrefix(out, left).tab();
		return cellPrefix(out, right);
	}

	private static void ratio(final ByteBuilder out, final long numerator, final long denominator){
		if(denominator==0){out.append("NA");}
		else{out.append(numerator/(double)denominator, 9);}
	}

	private static double seconds(final Cell cell, final int metric){
		assert(metric>=0 && metric<4) : "Timing reports have user, system, per-replicate CPU and elapsed seconds";
		return metric==0 ? cell.user : metric==1 ? cell.system : metric==2 ? cell.user+cell.system : cell.wall;
	}

	/** Input is a nonempty sorted observation array; CPU is summed before calling this. */
	private static double median(final double[] sorted){
		assert(sorted.length>0) : "The positive declared repeat count supplies median observations";
		final int middle=sorted.length/2;
		return (sorted.length&1)==1 ? sorted[middle] : sorted[middle-1]/2+sorted[middle]/2;
	}

	private void profileRow(final Output out, final byte[] name, final String kind, final long nanoseconds, final double processCpu){
		assert(nanoseconds>=0 && processCpu>0) : "Profile rows retain nonnegative worker clocks and their positive process denominator";
		final double seconds=nanoseconds/1e9;
		buffer.clear().append(name).tab().append(kind).tab().append(nanoseconds).tab().append(seconds, 9)
			.tab().append(processCpu, 6).tab().append(100*seconds/processCpu, 6)
			.append("\tCLOCK_OVERHEAD_INCLUDED_SHARED_PATH_AND_LOADING_NOT_ALLOCATED\n");
		out.writer.print(buffer);
	}

	private void check(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message+" at "+currentPath+":"+currentLine);}
	}

	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	private static boolean finite(final double value){return !Double.isNaN(value) && !Double.isInfinite(value);}
	private static byte[] ascii(final String value){return value.getBytes(StandardCharsets.US_ASCII);}

	/*--------------------------------------------------------------*/
	/*----------------       Owned observations     ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable lookup probe; only immutable copies enter the primitive-valued map. */
	private static final class Pool{
		int find(final LineParser1 parser, final int field){
			probe.set(parser, field);
			return index.get(probe);
		}
		int intern(final LineParser1 parser, final int field){
			probe.set(parser, field);
			return intern(probe);
		}
		int intern(final Target value){
			final int old=index.get(value);
			if(old>=0){return old;}
			require(value.to>value.from, "Native identity fields must be nonempty");
			final Target retained=value.copy();
			final int id=keys.size();
			keys.add(retained);
			index.put(retained, id);
			assert(index.size()==keys.size()) : "Each retained exact byte key owns one dense observation ID";
			return id;
		}
		int size(){return keys.size();}
		Target key(final int id){return keys.get(id);}
		final ObjectIntMap<Target> index=new ObjectIntMap<Target>(Target.class);
		final ArrayList<Target> keys=new ArrayList<Target>();
		final Target probe=new Target();
	}

	private static final class Locus{
		Locus(final int input_, final int family_){
			input=input_;
			family=family_;
		}
		final int input, family;
		int status=-1;
	}

	private static final class Grade{
		final long[] expected=new long[GRADE_COUNTS], observed=new long[GRADE_COUNTS];
	}

	/** Reused for every policy/repeat/worker comparison; no per-input boxed counters. */
	private static final class ComparisonScratch{
		ComparisonScratch(final int inputs, final int families, final int lineages){
			before=new int[inputs];
			after=new int[inputs];
			lost=new long[families];
			gained=new long[families];
			lineageLost=new long[Math.multiplyExact(families, lineages)];
			lineageGained=new long[lineageLost.length];
			lineageCounts=new long[Math.multiplyExact(lineages, LINEAGE_COUNTS)];
		}
		final int[] before, after;
		final long[] lost, gained;
		final long[] lineageLost, lineageGained, lineageCounts;
	}

	/** One retained result per declared invocation; no state is shared with callers. */
	private static final class Cell{
		Cell(final int policy_, final int worker_, final int repeat_){
			policy=policy_;
			worker=worker_;
			repeat=repeat_;
		}
		final int policy, worker, repeat;
		boolean timed;
		double user, system, wall;
		final IntHashMap inputRows=new IntHashMap();
		final IntList inputIds=new IntList(), references=new IntList(), tps=new IntList(), queryTps=new IntList();
		final LongList expectedQueries=new LongList();
		long[] actualQueries;
		boolean[] sourceSeen;
		final long[] sourceTotals=new long[4];
		final ArrayList<Grade> grades=new ArrayList<Grade>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String timings, inputManifest, summaries, raw, originalCounts, taxonomy, profile, profileTime, prefix;
	private String subtypeRaw, subtypeGrades;
	private NcrnaSubtypeComparison subtypeComparison;
	/** Optional declared class per selected input; labels never determine scoring. */
	private int[] lineageByInput, lineagePopulation;
	private long[] corpusInputs;
	private String baseline="shipped";
	private String[] policyNames;
	private int[] workerCounts={1, 16};
	private int repeats=3;
	/** Explicit mode: emits no grade/comparison files and makes no stability claim. */
	private boolean timingOnly;
	private final byte[][] policyBytes;
	private final int baselineIndex;
	private final Cell[] cells;
	private final PrintStream outstream;
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder(), keyBuffer=new ByteBuilder();
	private final Target probe=new Target();
	private final Pool inputs=new Pool(), families=new Pool(), parents=new Pool(), lineages=new Pool();
	private final ArrayList<Target> fastaPaths=new ArrayList<Target>(), referencePaths=new ArrayList<Target>();
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	private String currentPath="arguments";
	private long currentLine;
	private double overlap=-1;
	/** Counts differing cell comparisons once, even when several families differ. */
	private int unstableRepeats, unstableWorkers;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int INPUTS=0, TP=1, FP=2, FN=3, SCORED=4, MASKED_REF=5, MASKED_CALL=6,
		DUPLICATES=7, EXCLUDED=8, EXACT5=9, EXACT3=10, BOTH=11, GRADE_COUNTS=12;
	private static final String CONTEXT="policy\tthreads\trepeat\t";
	private static final String TAXONOMY_HEADER="input_id\tfasta\tdomain\tphylum\tcorpus_inputs\tcategory\ttaxonomy_status";
	private static final int[] TAXONOMY_LABEL_FIELDS={2, 3, 5, 6};
	private static final int LINEAGE_COUNTS=6;
	private static final String TIMING_HEADER="policy\tthreads\trepeat\tuser_seconds\tsystem_seconds\twall_seconds";
	private static final String ORIGINAL_HEADER="input_id\tsource_features\tCDS\ttRNA\trRNA\tprojected_RNA\tquery_features\tqualification";
	private static final String COMPARISON_HEADER="comparison\tleft_policy\tleft_workers\tleft_repeat\tright_policy\tright_workers\tright_repeat";
	private static final String PROFILE_HEADER="family\tstage\tinvocations\tcpu_ns\talignments\taligned_bases\thbm_calls\thbm_bases\tseed_hits\twindows";
	private static final String[] GLOBAL_STAGES={"SHARED_SEED_SCAN", "SSU_APPLY", "SSU_FINISH", "FIND_PATH", "TRNA"};
}
