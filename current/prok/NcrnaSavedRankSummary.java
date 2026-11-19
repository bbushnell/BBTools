package prok;

import java.nio.charset.StandardCharsets;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Streams the saved locus minima and parent inventory; never opens sequences.
 * Derived from template.A_SampleByteFile's parser/native byte-I/O lifecycle.
 * Sums are over common loci, not measured alignments from a literal caller walk.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSavedRankSummary{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSavedRankSummary.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			int methods=148, models=508, loci=2156, common=2007, zero=3, gaps=110;
			int[] selected={0, 54, 57, 84};
			String strata=null, panels=null;
			int groups=0;
			boolean fullCapture=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("methods")){selected=Parse.parseIntArray(b, ",");}
				else if(a.equals("strata")){strata=b;}
				else if(a.equals("panels")){panels=b;}
				else if(a.equals("expectedgroups")){groups=Parse.parseIntKMG(b);}
				else if(a.equals("fullcapture")){fullCapture=Parse.parseBoolean(b);}
				else if(a.equals("expectedmethods")){methods=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmodels")){models=Parse.parseIntKMG(b);}
				else if(a.equals("expectedloci")){loci=Parse.parseIntKMG(b);}
				else if(a.equals("expectedcommon")){common=Parse.parseIntKMG(b);}
				else if(a.equals("expectedzero")){zero=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgaps")){gaps=Parse.parseIntKMG(b);}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown saved-rank parameter: "+arg);}
			}
			require((parser.in1!=null ^ panels!=null) && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use exactly one of in=saved_rank_prefix or panels=manifest.tsv, with out=fresh_prefix and no overwrite, append or sampling");
			require(!fullCapture || panels!=null, "fullcapture requires explicit complete source panels, never a single unbound rank prefix");
			if(panels==null){process(parser.in1, parser.out1, strata, selected, methods, models, loci, common, zero, gaps);}
			else{
				require(strata==null, "Combined panels declare their own complete stratum-label paths");
				combine(panels, parser.out1, selected, methods, models, groups, fullCapture);
			}
			pp.outstream.println("NCRNA_SAVED_RANK_SUMMARY_COMPLETE_NO_CALLER_OR_ALIGNMENT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Expected counts bind this report to the reviewed complete saved experiment. */
	public static void process(final String input, final String output, final String strata, final int[] selected,
			final int methods, final int models, final int loci, final int common, final int zero, final int gaps){
		require(input!=null && output!=null && selected!=null && selected.length>1 && selected[0]==0
			&& methods>1 && models>0 && loci>0 && common>0 && common<=loci && zero>=0 && gaps>=0,
			"Saved rank reporting requires a complete locus inventory, positive common set and native method0 first");
		final String rankPath=input+".loci.tsv.gz", parentPath=input+".parents.tsv.gz";
		final String statsPath=output+".summary.tsv", zeroPath=output+".zero_support.tsv", gapPath=output+".no_compatible.tsv";
		require(Tools.testInputFiles(false, true, rankPath, parentPath)
			&& (strata==null || Tools.testInputFiles(false, true, strata))
			&& Tools.testOutputFiles(false, false, false, statsPath, zeroPath, gapPath)
			&& Tools.testForDuplicateFiles(true, rankPath, parentPath, strata, statsPath, zeroPath, gapPath),
			"Saved input paths and fresh outputs must be readable, distinct and nonoverwriting");
		new NcrnaSavedRankSummary(selected, methods, models).run(rankPath, parentPath, strata,
			statsPath, zeroPath, gapPath, loci, common, zero, gaps);
	}

	private NcrnaSavedRankSummary(final int[] selected_, final int methods_, final int models_){
		require(selected_!=null && selected_.length>1 && selected_[0]==0 && methods_>1 && models_>0,
			"Saved ranking requires native method0 first and complete positive model/method domains");
		int previous=-1;
		for(final int method:selected_){
			require(method>previous && method<methods_, "Selected method IDs must be distinct, increasing and present in the saved grid");
			previous=method;
		}
		selected=selected_.clone();
		methods=methods_;
		models=models_;
		ranks=new int[selected.length];
		sums=new long[STRATA.length][selected.length];
		histograms=new long[STRATA.length][selected.length][models+1];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Saved Rows          ----------------*/
	/*--------------------------------------------------------------*/

	/** Combine disjoint assembly groups, retaining all source validation and no pooled stratum. */
	public static void combine(final String panels, final String output, final int[] selected,
			final int methods, final int models, final int expectedGroups){
		combine(panels, output, selected, methods, models, expectedGroups, false);
	}

	/** Full capture keeps contributors and the first captured representative of aliases.
	 * A zero expectedGroups requests a census of the complete, explicit source panels. */
	public static void combine(final String panels, final String output, final int[] selected,
			final int methods, final int models, final int expectedGroups, final boolean fullCapture){
		require(panels!=null && output!=null && (expectedGroups>0 || fullCapture && expectedGroups==0),
			"Combined reporting requires a panel manifest and the complete selected assembly-group count");
		RrnaResourceIO.fresh(output, ".summary.tsv", ".zero_support.tsv", ".no_compatible.tsv", ".inputs.tsv");
		if(fullCapture){RrnaResourceIO.fresh(output, ".panels.tsv", ".aliases.tsv", ".population.tsv");}
		final NcrnaSavedRankSummary tool=new NcrnaSavedRankSummary(selected, methods, models);
		tool.fullCapture=fullCapture;
		tool.combined(panels, output, expectedGroups);
	}

	private void combined(final String manifest, final String output, final int expectedGroups){
		final LineParser1 table=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		try(Input source=new Input(manifest); Output stats=new Output(output+".summary.tsv");
			Output zero=new Output(output+".zero_support.tsv"); Output gaps=new Output(output+".no_compatible.tsv");
			Output inventory=new Output(output+".inputs.tsv");
			Output effective=fullCapture ? new Output(output+".panels.tsv") : null;
			Output aliases=fullCapture ? new Output(output+".aliases.tsv") : null;
			Output population=fullCapture ? new Output(output+".population.tsv") : null){
			header(source, PANEL_HEADER);
			zero.writer.println("panel\t"+PARENT_HEADER+"\tnative_rank\tselected_methods\tselected_ranks\tinterpretation");
			gaps.writer.println("panel\t"+PARENT_HEADER+"\tinterpretation");
			inventory.writer.println("panel\tinput_slot\tinput_id\tfasta\tresolved_assembly_group\ttarget_loci\tcommon_loci\tnative_zero_support\tno_compatible_model\tno_window\tambiguous_only");
			if(fullCapture){
				effective.writer.println(PANEL_HEADER);
				aliases.writer.println("panel\tinput_slot\tinput_id\tfasta\tassembly_group\tdisposition");
				population.writer.println("captured_inputs\tdistinct_assembly_groups\tdiscarded_aliases");
			}
			int candidates=0, duplicates=0;
			for(byte[] row=source.reader.nextLine(); row!=null; row=source.reader.nextLine()){
				table.set(row);
				require(table.terms()==11, "Every panel row must bind source paths, selection and complete source counts");
				panelName=text(table, 0);
				require(panelName.matches("[A-Za-z0-9_.-]+") && names.add(panelName), "Each panel needs a distinct printable name");
				final String prefix=text(table, 1), strata=text(table, 2), assemblies=text(table, 3), include=text(table, 4);
				require(prefix.startsWith("/") && strata.startsWith("/") && assemblies.startsWith("/")
					&& (include.equals(".") || include.startsWith("/")), "Panel paths must be explicit absolute saved inputs; dot selects every assembly row");
				final int inputCount=table.parseInt(5), selectedCount=table.parseInt(6);
				require(!fullCapture || include.equals(".") && selectedCount==inputCount,
					"Full capture requires every original input as a candidate; no contributor/subset filter is permitted");
				final int loci=table.parseInt(7), common=table.parseInt(8), zeroCount=table.parseInt(9), gapCount=table.parseInt(10);
				bindings(prefix+".inputs.tsv", assemblies, include, inputCount, selectedCount);
				candidates+=inputCount;
				int kept=0;
				for(final BoundInput input:inputs){if(input.selected){kept++;}}
				if(fullCapture && kept>0){
					buffer.clear();
					for(int column=0; column<11; column++){
						if(column>0){buffer.tab();}
						if(column==6){buffer.append(kept);}
						else{table.setBounds(column); buffer.append(table.line(), table.a(), table.b()-table.a());}
					}
					effective.writer.print(buffer.nl());
				}
				scan(prefix+".loci.tsv.gz", prefix+".parents.tsv.gz", strata, zero.writer, gaps.writer, loci, common, zeroCount, gapCount);
				for(int slot=0; slot<inputs.length; slot++){
					final BoundInput input=inputs[slot];
					if(!input.selected){
						if(fullCapture){
							aliases.writer.print(buffer.clear().append(panelName).tab().append(slot).tab().append(input.idBytes)
								.tab().append(input.pathBytes).tab().append(input.group).append("\tDUPLICATE_ASSEMBLY_FIRST_CAPTURE_KEPT\n"));
							duplicates++;
						}
						continue;
					}
					assert(input.targets==input.common+input.zero+input.gaps+input.noWindow+input.ambiguous)
						: "Selected input inventory partitions every target using the saved producer dispositions";
					inventory.writer.print(buffer.clear().append(panelName).tab().append(slot).tab().append(input.idBytes)
						.tab().append(input.pathBytes).tab().append(input.group).tab().append(input.targets).tab().append(input.common)
						.tab().append(input.zero).tab().append(input.gaps).tab().append(input.noWindow).tab().append(input.ambiguous).nl());
				}
			}
			require(names.size()>0 && selectedGroups.size()>0 && (expectedGroups==0 || selectedGroups.size()==expectedGroups),
				"Combined inputs must contain every declared disjoint assembly group, including zero-target groups; actual="+selectedGroups.size());
			assert(!fullCapture || candidates==selectedGroups.size()+duplicates)
				: "Full capture partitions every original input into one retained assembly representative or a declared alias";
			report(stats.writer, true);
			if(fullCapture){population.writer.print(buffer.clear().append(candidates).tab().append(selectedGroups.size()).tab().append(duplicates).nl());}
		}
	}

	/** Selection binds dense original slots, original IDs/paths and actual noncontribution. */
	private void bindings(final String captured, final String assemblies, final String include, final int inputCount, final int selectedCount){
		require(inputCount>0 && selectedCount>0 && selectedCount<=inputCount,
			"A panel must declare its complete source population and nonempty selected subset");
		inputs=new BoundInput[inputCount];
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), paths=new ObjectSet<String>(String.class);
		try(Input in=new Input(assemblies)){
			header(in, NcrnaDdlLocusOverlap.ASSEMBLY_HEADER);
			int count=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7 && count<inputs.length && p.parseInt(0)==count,
					"Assembly metadata must retain every dense original captured input slot");
				final String id=text(p, 1), path=text(p, 2);
				require(!id.isEmpty() && path.startsWith("/") && ids.add(id) && paths.add(path),
					"Assembly census IDs and paths must uniquely identify every original input");
				inputs[count++]=new BoundInput(id, path, p.parseInt(3), p.parseLong(4), p.termEquals("NO_KNOWN_MEMBER_CONTRIBUTOR", 5));
			}
			require(count==inputCount, "The assembly census cannot omit zero-target inputs");
		}
		try(Input in=new Input(captured)){
			header(in, NcrnaDdlTrainingOverlap.SAMPLE_HEADER);
			int count=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==13 && count<inputs.length && p.parseInt(0)==count
					&& same(p, 1, inputs[count].idBytes) && same(p, 2, inputs[count].pathBytes)
					&& p.parseInt(3)>0 && p.parseLong(4)>0,
					"Assembly rows must match every original capture slot/ID/path, including zero-target inputs");
				count++;
			}
			require(count==inputs.length, "The original capture accounting must retain the complete panel population");
		}
		int selected=0;
		if(include.equals(".")){
			for(final BoundInput input:inputs){
				select(input);
				selected++;
			}
		}else{
			try(Input in=new Input(include)){
				header(in, REUSE_HEADER);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					p.set(row);
					require(p.terms()==4, "Reused selection requires the original group, captured slot, ID and path");
					final int slot=p.parseInt(1);
					require(slot>=0 && slot<inputs.length, "A reused captured slot must be inside its original complete input census");
					final BoundInput input=inputs[slot];
					require(p.parseInt(0)==input.group && same(p, 2, input.idBytes) && same(p, 3, input.pathBytes),
						"Reused groups must retain the exact original captured slot, input ID and FASTA path");
					select(input);
					selected++;
				}
			}
		}
		require(selected==selectedCount, "The selected subset must contain every declared reused/new group exactly once; actual="+selected);
	}

	private void select(final BoundInput input){
		require(!input.selected && input.group>=0 && (fullCapture || input.votes==0 && input.noncontributor),
			"Only resolved noncontributing assembly groups can enter the expanded panel");
		if(fullCapture && selectedGroups.contains(input.group)){return;}
		require(!selectedGroups.contains(input.group) && selectedIds.add(input.id) && selectedPaths.add(input.path),
			"Reused and new panels must be disjoint in assembly groups, original input IDs and FASTA paths");
		selectedGroups.put(input.group, 1);
		input.selected=true;
	}

	private BoundInput boundInput(final LineParser1 p){
		final int slot=p.parseInt(0);
		require(slot>=0 && slot<inputs.length, "Every saved parent must belong to its original complete input census");
		final BoundInput input=inputs[slot];
		require(same(p, 1, input.idBytes) && same(p, 2, input.pathBytes),
			"Parent input ID/path must agree with the original assembly census before subset selection");
		return input;
	}

	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 cannot change panel or assembly identity");
		return value;
	}

	private static boolean same(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			if(p.line()[p.a()+i]!=expected[i]){return false;}
		}
		return true;
	}

	private void run(final String rankPath, final String parentPath, final String strataPath, final String statsPath,
			final String zeroPath, final String gapPath, final int expectedLoci, final int expectedCommon,
			final int expectedZero, final int expectedGaps){
		try(Output stats=new Output(statsPath); Output zero=new Output(zeroPath); Output gaps=new Output(gapPath)){
			zero.writer.println(PARENT_HEADER+"\tnative_rank\tselected_methods\tselected_ranks\tinterpretation");
			gaps.writer.println(PARENT_HEADER+"\tinterpretation");
			scan(rankPath, parentPath, strataPath, zero.writer, gaps.writer, expectedLoci, expectedCommon, expectedZero, expectedGaps);
			report(stats.writer, strataPath!=null);
		}
	}

	/** Validate every source row, even when only a declared subset contributes. */
	private void scan(final String rankPath, final String parentPath, final String strataPath,
			final ByteStreamWriter zero, final ByteStreamWriter gaps, final int expectedLoci,
			final int expectedCommon, final int expectedZero, final int expectedGaps){
		require(expectedLoci>=0 && expectedCommon>=0 && expectedZero>=0 && expectedGaps>=0,
			"Each panel must declare its complete source counts before subset selection");
		all=fullCommon=zeroCount=gapCount=0;
		seen.clear();
		try(Input source=new Input(rankPath); Input parents=new Input(parentPath);
			Input strata=strataPath==null ? null : new Input(strataPath)){
			header(source, LOCUS_HEADER);
			header(parents, PARENT_HEADER);
			if(strata!=null){header(strata, STRATA_HEADER);}
			for(byte[] first=source.reader.nextLine(); first!=null; first=source.reader.nextLine()){
				parentRow=nextTarget(parents);
				require(parentRow!=null, "Saved locus rows exceed the original target-parent inventory");
				parent.set(parentRow);
				inputSlot=parent.parseInt(0);
				inputId=parent.parseString(1);
				parentId=parent.parseString(3);
				currentInput=inputs==null ? null : boundInput(parent);
				stratum=readStratum(strata);
				require(seen.add(parentId), "A saved target parent must occur in exactly one complete method block: "+parentId);
				int picked=0;
				for(int method=0; method<methods; method++){
					final byte[] row=method==0 ? first : source.reader.nextLine();
					require(row!=null, "Truncated method block for target parent "+parentId);
					line.set(row);
					require(line.terms()==12 && line.parseInt(0)==inputSlot && line.termEquals(inputId, 1)
						&& line.termEquals(parentId, 2) && line.parseInt(3)==method,
						"Saved locus rows must preserve parent identity and every method in grid order: "+parentId);
					final int rank=line.parseInt(8);
					require(rank==-1 || rank>0 && rank<=models, "Saved rank is outside the declared model order: "+parentId);
					disposition(method, rank);
					if(picked<selected.length && method==selected[picked]){ranks[picked++]=rank;}
				}
				require(picked==selected.length, "Each selected method must be present in the completed locus block");
				consume(zero, gaps);
				all++;
			}
			require(nextTarget(parents)==null && (strata==null || strata.reader.nextLine()==null)
				&& all==expectedLoci && fullCommon==expectedCommon
				&& zeroCount==expectedZero && gapCount==expectedGaps,
				"Saved report must conserve reviewed target/common/zero-support/library-gap counts: actual="
				+all+","+fullCommon+","+zeroCount+","+gapCount);
		}
	}

	/** Mirrors RankNcrnaWindows.reason; the four window counters are method-invariant. */
	private void disposition(final int method, final int rank){
		final int windows=line.parseInt(4), compatible=line.parseInt(5), ambiguous=line.parseInt(6), unique=line.parseInt(7);
		require(windows>=0 && compatible>=0 && compatible<=windows && ambiguous>=0 && unique>=0
			&& (long)ambiguous+unique==compatible,
			"Saved window dispositions must partition compatible events into unique and ambiguous: "+parentId);
		if(method==0){
			windowCounts[0]=windows;
			windowCounts[1]=compatible;
			windowCounts[2]=ambiguous;
			windowCounts[3]=unique;
		}else{
			require(windows==windowCounts[0] && compatible==windowCounts[1] && ambiguous==windowCounts[2] && unique==windowCounts[3],
				"A ranking method cannot change the saved window/oracle dispositions: "+parentId);
		}
		final String expected=windows==0 ? "NO_CAPTURED_WINDOW" : compatible==0 ? "NO_COMPATIBLE_MODEL"
			: unique==0 ? "AMBIGUOUS_SPANS_ONLY" : method==0 && rank<0 ? "ZERO_NATIVE_INDEX_SUPPORT" : "RANKED";
		require(line.termEquals(expected, 9) && (expected.equals("RANKED") ? rank>0 : rank==-1),
			"Saved rank and reason must agree with the producer's window/oracle state for every method: "+parentId+", method="+method);
		if(method==0){reason=expected;}
	}

	private void consume(final ByteStreamWriter zero, final ByteStreamWriter gaps){
		if(ranks[0]>0){fullCommon++;}
		else if(reason.equals("ZERO_NATIVE_INDEX_SUPPORT")){zeroCount++;}
		else if(reason.equals("NO_COMPATIBLE_MODEL")){gapCount++;}
		// Validation includes excluded rows; only the selected inventory enters sums.
		for(int i=1; i<ranks.length; i++){
			if(ranks[0]>0 || reason.equals("ZERO_NATIVE_INDEX_SUPPORT")){
				require(ranks[i]>0, "Every unique-oracle locus must retain the requested sketch rank, including excluded inputs");
			}
		}
		if(currentInput!=null){
			if(!currentInput.selected){return;}
			currentInput.targets++;
			if(ranks[0]>0){currentInput.common++;}
			else if(reason.equals("ZERO_NATIVE_INDEX_SUPPORT")){currentInput.zero++;}
			else if(reason.equals("NO_COMPATIBLE_MODEL")){currentInput.gaps++;}
			else if(reason.equals("NO_CAPTURED_WINDOW")){currentInput.noWindow++;}
			else if(reason.equals("AMBIGUOUS_SPANS_ONLY")){currentInput.ambiguous++;}
		}
		if(ranks[0]>0){
			require(reason.equals("RANKED"), "A positive native rank must retain its saved RANKED disposition");
			for(int i=0; i<ranks.length; i++){
				require(ranks[i]>0, "A sketch must rank every native-rankable unique-oracle locus before a common-set mean is defined");
				sums[0][i]+=ranks[i];
				sums[stratum][i]+=ranks[i];
				histograms[0][i][ranks[i]]++;
				histograms[stratum][i][ranks[i]]++;
			}
			common[0]++;
			common[stratum]++;
		}else if(reason.equals("ZERO_NATIVE_INDEX_SUPPORT")){
			beginParent().tab().append(ranks[0]).tab();
			for(int i=1; i<selected.length; i++){
				if(i>1){buffer.append(',');}
				buffer.append(selected[i]);
			}
			buffer.tab();
			for(int i=1; i<ranks.length; i++){
				require(ranks[i]>0, "Zero-native-support cases must retain the requested sketch ranks");
				if(i>1){buffer.append(',');}
				buffer.append(ranks[i]);
			}
			buffer.append("\tORACLE_ACCEPTABLE_WITH_ZERO_NATIVE_INDEX_SUPPORT_NOT_DETECTION_RECOVERY\n");
			zero.print(buffer);
		}else if(reason.equals("NO_COMPATIBLE_MODEL")){
			beginParent().append("\tNO_COMPATIBLE_ACCEPTED_MODEL_IN_CAPTURED_WINDOWS;CAUSE_UNRESOLVED_WINDOW_VERIFIER_OR_LIBRARY\n");
			gaps.print(buffer);
		}else{
			require(reason.equals("NO_CAPTURED_WINDOW") || reason.equals("AMBIGUOUS_SPANS_ONLY"),
				"Every unranked locus needs a recognized saved disposition: "+reason);
		}
	}

	private ByteBuilder beginParent(){
		assert(parentRow!=null) : "Named unresolved cases retain the bound original parent row";
		buffer.clear();
		if(panelName!=null){buffer.append(panelName).tab();}
		return buffer.append(parentRow);
	}

	/** Unknown is explicit; absence of a label file never means independent holdout. */
	private int readStratum(final Input in){
		if(in==null){return 3;}
		final byte[] row=in.reader.nextLine();
		require(row!=null, "Declared seen-status labels must cover every target parent");
		labels.set(row);
		require(labels.terms()==4 && labels.parseInt(0)==inputSlot && labels.termEquals(inputId, 1)
			&& labels.termEquals(parentId, 2), "Seen-status labels must preserve the complete saved parent order and identities");
		for(int s=1; s<STRATA.length; s++){
			if(labels.termEquals(STRATA[s], 3)){return s;}
		}
		throw new IllegalArgumentException("Seen-status labels must be SEEN, NOT_SEEN or UNKNOWN: "+parentId);
	}

	/** Reuses the original report's nearest-rank quantile definition. */
	private void report(final ByteStreamWriter out, final boolean labelled){
		out.println("method\tcommon_loci\tsum_of_locus_minimum_ranks\tmean_locus_minimum_rank\tinterpretation\tseen_stratum\twithin10\twithin50\twithin100\tp90\tp95\tp99");
		final int[] order=panelName==null ? new int[]{0, 1, 2, 3} : new int[]{2, 1, 3};
		for(final int s:order){
			for(int i=0; i<selected.length; i++){
				final long[] histogram=histograms[s][i];
				buffer.clear().append(selected[i]).tab().append(common[s]).tab().append(sums[s][i]).tab();
				if(common[s]==0){buffer.append("NA");}else{buffer.append((double)sums[s][i]/common[s], 9);}
				buffer.append("\tSAVED_LOCUS_MINIMUM_PROXY_NOT_MEASURED_CALLER_ALIGNMENTS;")
					.append(labelled ? "DECLARED_INPUT_STRATA" : "SEEN_STATUS_NOT_ASSESSED")
					.tab().append(STRATA[s]).tab().append(covered(histogram, 10)).tab()
					.append(covered(histogram, 50)).tab().append(covered(histogram, 100)).tab()
					.append(RankNcrnaWindows.quantile(histogram, 0.90)).tab()
					.append(RankNcrnaWindows.quantile(histogram, 0.95)).tab()
					.append(RankNcrnaWindows.quantile(histogram, 0.99)).nl();
				out.print(buffer);
			}
		}
	}

	private static long covered(final long[] histogram, final int limit){
		assert(histogram!=null && histogram.length>1 && limit>0) : "Coverage shares the original positive-rank histogram domain";
		long sum=0;
		for(int i=1; i<=limit && i<histogram.length; i++){sum+=histogram[i];}
		return sum;
	}

	/** Parent inventory contains other family/mask rows; preserve its target order. */
	private byte[] nextTarget(final Input in){
		for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
			parent.set(row);
			require(parent.terms()==14, "Saved parent inventory must retain its complete source/coordinate columns");
			if(inputs!=null){boundInput(parent);}
			if(parent.termEquals("true", 9)){return row;}
			require(parent.termEquals("false", 9), "Saved target flag must be an explicit boolean");
		}
		return null;
	}

	private static void header(final Input in, final String expected){
		final byte[] row=in.reader.nextLine();
		require(row!=null && new String(row, java.nio.charset.StandardCharsets.UTF_8).equals(expected),
			"Saved rank schemas must match the reviewed RankNcrnaWindows report before fields are interpreted");
	}

	private static void require(final boolean condition, final String why){
		if(!condition){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] selected, ranks;
	private final int[] windowCounts=new int[4];
	private final long[][] sums;
	private final long[][][] histograms;
	private final int[] common=new int[STRATA.length];
	private final int methods, models;
	private final LineParser1 line=new LineParser1('\t'), parent=new LineParser1('\t');
	private final LineParser1 labels=new LineParser1('\t');
	private final ObjectSet<String> seen=new ObjectSet<String>(String.class);
	private final ByteBuilder buffer=new ByteBuilder();
	private byte[] parentRow;
	private String inputId, parentId, reason;
	private int inputSlot, all, fullCommon, zeroCount, gapCount, stratum;
	private BoundInput[] inputs;
	private BoundInput currentInput;
	private String panelName;
	private boolean fullCapture;
	private final IntHashMap selectedGroups=new IntHashMap();
	private final ObjectSet<String> selectedIds=new ObjectSet<String>(String.class), selectedPaths=new ObjectSet<String>(String.class);
	private static final class BoundInput{
		BoundInput(final String id_, final String path_, final int group_, final long votes_, final boolean noncontributor_){
			require(group_>=-1 && votes_>=0, "Assembly metadata preserves resolved/unknown group and nonnegative kept-vote count");
			id=id_;
			path=path_;
			idBytes=id.getBytes(StandardCharsets.UTF_8);
			pathBytes=path.getBytes(StandardCharsets.UTF_8);
			group=group_;
			votes=votes_;
			noncontributor=noncontributor_;
		}
		final String id, path;
		final byte[] idBytes, pathBytes;
		final int group;
		final long votes;
		final boolean noncontributor;
		boolean selected;
		int targets, common, zero, gaps, noWindow, ambiguous;
	}
	private static final String[] STRATA={"ALL_COMMON", "SEEN", "NOT_SEEN", "UNKNOWN"};
	public static final String STRATA_HEADER="input_slot\tinput_id\tparent\tseen_stratum";
	public static final String LOCUS_HEADER="input_slot\tinput_id\tparent\tmethod\twindows\tcompatible_windows\tambiguous_only_windows\tunique_windows\tbest_rank\treason\tany_unique_model_native_shortlist_eligible\tany_uniquely_attributed_oracle_model_actually_attempted";
	public static final String PARENT_HEADER="input_slot\tinput_id\tfasta\tparent\tseqid\tstart1\tstop1\tstrand\tmask\ttarget\tprior\tmulti_supported\tscored_models\tmarginal_models";
	public static final String PANEL_HEADER="panel\trank_prefix\tstrata\tassemblies\tinclude\texpected_inputs\texpected_selected_inputs\texpected_loci\texpected_common\texpected_zero\texpected_gaps";
	public static final String REUSE_HEADER="resolved_assembly_group\tcaptured_slot\tinput_id\tfasta";
}
