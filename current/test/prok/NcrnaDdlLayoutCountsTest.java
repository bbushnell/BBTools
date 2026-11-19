package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser1;
import prok.NcrnaDdlLayoutCounts;
import prok.NcrnaDdlSavedData;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.OptimizeNcrnaPartitionWeights;
import structures.ByteBuilder;

/** Synthetic interchange tests reuse the saved two-model fixture. Its metadata
 * describes literal toy arrays, never biological resource or memory measurements.
 * @author Ganyu
 */
final class NcrnaDdlLayoutCountsTest{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	static void run(final String root, final OptimizeNcrnaPartitionWeights.Config original){
		final OptimizeNcrnaPartitionWeights.Config c=config(original, root+"/layouts_fit");
		c.layouts=root+"/layouts.tsv";
		c.counts=root+"/layout_counts.tsv.gz";
		metadata(c.layouts);
		counts(original.output+".counts.tsv.gz", c.counts);
		OptimizeNcrnaPartitionWeights.run(c);
		final NcrnaDdlLayoutCounts historical=new NcrnaDdlLayoutCounts(c.layouts, 2);
		check(historical.historical && historical.layouts[0].terms==3 && historical.layouts[1].terms==4,
			"Version1 remains explicitly historical and preserves its three/four-term interpretation");
		selection(root, original, c);
		try(Input a=new Input(c.layouts); Input b=new Input(c.output+".layouts.tsv")){
			for(byte[] row; (row=a.reader.nextLine())!=null;){
				check(Arrays.equals(row, b.reader.nextLine()), "Validated layout payload metadata must survive byte-for-byte");
			}
			check(b.reader.nextLine()==null, "Metadata copy cannot add rows");
		}
		for(final String fault:new String[]{"missing", "extra", "identity", "overflow", "forbidden"}){
			final String path=root+"/layout_"+fault+".tsv";
			faultCounts(c.counts, path, fault);
			final OptimizeNcrnaPartitionWeights.Config bad=config(original, root+"/reject_layout_"+fault);
			bad.layouts=c.layouts;
			bad.counts=path;
			bad.fitLayouts=c.control;
			boolean rejected=false;
			try{OptimizeNcrnaPartitionWeights.run(bad);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected && !new File(bad.output+".trials.tsv.gz").exists(),
				"Malformed layout counts must fail before publishing trials: "+fault);
		}
		final String duplicate=root+"/duplicate_layouts.tsv";
		try(Input input=new Input(c.layouts); Output output=new Output(duplicate)){
			byte[] last=null;
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				output.writer.println(row);
				last=row;
			}
			output.writer.println(last);
		}
		boolean rejected=false;
		try{new NcrnaDdlLayoutCounts(duplicate, 2);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Duplicate layout identifiers must not overwrite or merge independent fits");
		for(final boolean inKey:new boolean[]{false, true}){
			final String path=root+"/zero_member_origins_"+inKey+".tsv";
			try(Input input=new Input(c.layouts); Output output=new Output(path)){
				output.writer.println(input.reader.nextLine());
				final String[] fields=new String(input.reader.nextLine(), StandardCharsets.UTF_8).split("\t");
				fields[14]="0";
				if(inKey){fields[10]="LOW_KEY_BIT";}
				output.writer.println(String.join("\t", fields).getBytes(StandardCharsets.UTF_8));
				for(byte[] row; (row=input.reader.nextLine())!=null;){output.writer.println(row);}
			}
			rejected=false;
			try{new NcrnaDdlLayoutCounts(path, 2);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected!=inKey,
				"An in-key discriminator permits zero stored origins; a mixed sidecar cannot hide required origin bytes");
		}
		current(root, original);
		System.err.println("NCRNA_DDL_LAYOUT_COUNTS_TEST_PASS");
	}

	/** Selecting one completed scientific search must preserve every result byte
	 * while leaving unselected fits absent. Validation still covers ALL counts. */
	private static void selection(final String root, final OptimizeNcrnaPartitionWeights.Config original,
			final OptimizeNcrnaPartitionWeights.Config all){
		for(final String id:IDS){
			final OptimizeNcrnaPartitionWeights.Config selected=config(original, root+"/selected_"+id);
			selected.layouts=all.layouts;
			selected.counts=all.counts;
			selected.fitLayouts=id;
			if(id.equals(selected.control)){
				OptimizeNcrnaPartitionWeights.main(new String[]{"layouts="+selected.layouts, "counts="+selected.counts,
					"control="+selected.control, "fitlayouts="+id, "panels="+selected.panels,
					"selection="+selected.selection, "ref="+selected.references, "identities="+selected.identities,
					"ranks="+selected.ranks, "out="+selected.output, "expectedmodels="+selected.models,
					"expectedgroups="+selected.groups, "randomdraws="+selected.randomDraws, "t=1"});
			}else{OptimizeNcrnaPartitionWeights.run(selected);}
			for(final String suffix:FIT_SUFFIXES){
				compareSelected(all.output+suffix, selected.output+suffix, id);
				final LineParser1 p=new LineParser1('\t');
				try(Input input=new Input(selected.output+suffix)){
					input.reader.nextLine();
					for(byte[] row; (row=input.reader.nextLine())!=null;){
						p.set(row);
						check(p.termEquals(id, 0), "Unselected layouts must never be fitted or reported");
					}
				}
			}
			timing(selected.output+".timing.tsv", new String[]{id});
		}
		final OptimizeNcrnaPartitionWeights.Config remaining=config(original, root+"/selected_remaining");
		remaining.layouts=all.layouts;
		remaining.counts=all.counts;
		remaining.fitLayouts="mixed256,long256";
		OptimizeNcrnaPartitionWeights.run(remaining);
		timing(remaining.output+".timing.tsv", new String[]{"long256", "mixed256"});
		int ordinal=0;
		for(final String ids:new String[]{"absent", "long256,long256", "", ",long256", "long256,"}){
			final OptimizeNcrnaPartitionWeights.Config invalid=config(original, root+"/selection_reject_"+ordinal++);
			invalid.layouts=all.layouts;
			invalid.counts=all.counts;
			invalid.fitLayouts=ids;
			boolean rejected=false;
			try{OptimizeNcrnaPartitionWeights.run(invalid);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected && !new File(invalid.output+".trials.tsv.gz").exists(),
				"Unknown, duplicate or empty selection tokens must fail before publishing trials");
		}
		System.err.println("NCRNA_DDL_SELECTED_LAYOUT_FIT_TEST_PASS");
	}

	private static void timing(final String path, final String[] ids){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			check(input.reader.nextLine()!=null, "Per-layout timing requires its explicit units/scope header");
			for(final String id:ids){
				final byte[] row=input.reader.nextLine();
				check(row!=null, "Each selected layout must have exactly one completed timing row");
				p.set(row);
				final int terms=id.equals("long256") ? 3 : 4;
				check(p.terms()==7 && p.termEquals(id, 0) && p.parseInt(2)==terms
					&& p.parseInt(3)>0 && p.parseInt(3)<=(terms==3 ? 471 : 1011)
					&& Double.isFinite(p.parseDouble(4)) && p.parseDouble(4)>=0
					&& Double.isFinite(p.parseDouble(5)) && p.parseDouble(5)>=0
					&& p.termEquals("FIT_AND_REPORT_ENQUEUE_EXCLUDES_INPUT_AND_CONTROL_VALIDATION", 6),
					"Timing must retain manifest order, actual bounded score work and finite current-thread/elapsed seconds");
			}
			check(input.reader.nextLine()==null, "Skipped layouts cannot emit timing rows or repeat the control fit");
		}
	}

	private static OptimizeNcrnaPartitionWeights.Config config(final OptimizeNcrnaPartitionWeights.Config original, final String output){
		final OptimizeNcrnaPartitionWeights.Config c=new OptimizeNcrnaPartitionWeights.Config();
		c.panels=original.panels;
		c.selection=original.selection;
		c.references=original.references;
		c.identities=original.identities;
		c.ranks=original.ranks;
		c.output=output;
		c.models=original.models;
		c.groups=original.groups;
		c.randomDraws=original.randomDraws;
		c.control="legacy_member384_consensus128_both";
		return c;
	}

	private static void metadata(final String path){
		try(Output output=new Output(path)){
			output.writer.println(NcrnaDdlLayoutCounts.LAYOUT_HEADER);
			for(final String id:IDS){
				final boolean legacy=id.equals(IDS[1]);
				final int member=legacy ? 384 : 256, consensus=legacy ? 128 : 256;
				final int slots=member+consensus;
				final ByteBuilder b=new ByteBuilder();
				b.append(id).append("\tTOY\t").append(legacy ? "PARTITIONED" : "MIXED").tab().append(member)
					.tab().append(legacy ? 128 : 0).tab().append(legacy ? 256 : 0)
					.tab().append(id.equals("long256") ? "LONG_ONLY" : "BOTH").tab().append(consensus)
					.append(legacy ? "\t4\t12\tLEGACY_PARTITION_ORIGINS\tWRAP" : "\t5\t10\tSIDECAR_ORIGIN\tCLAMP")
					.append("\t2").tab().append(4L*member).tab().append(legacy ? 0 : 2L*member)
					.tab().append(4L*consensus).tab().append(2L*consensus)
					.tab().append(slots).tab().append(slots).append("\t0\t4").tab().append(4L*(slots+4))
					.tab().append(8L*slots).tab().append(2L*slots).append("\t4\t0\t0\n");
				output.writer.print(b);
			}
		}
	}

	/** Deliberately put the control second: readers must use its explicit ID. */
	private static void counts(final String source, final String path){
		final LineParser1 p=new LineParser1('\t');
		try(Output output=new Output(path)){
			output.writer.println(NcrnaDdlLayoutCounts.COUNT_HEADER);
			for(final String id:IDS){
				final String geometry=id.equals("long256") ? "CONSENSUS10_ONLY" : "CONSENSUS_BOTH_LENGTHS";
				try(Input input=new Input(source)){
					input.reader.nextLine();
					for(byte[] row; (row=input.reader.nextLine())!=null;){
						p.set(row);
						check(p.terms()==14, "The original two-model capture requires six identity fields and four counts per model");
						if(!p.termEquals(geometry, 0) || p.parseInt(1)!=512){continue;}
						final ByteBuilder out=new ByteBuilder().append(id);
						identity(out, p);
						// Synthetic historical BOTH rows exercise an active fourth
						// term while preserving the unit total of the saved tiny fixture.
						for(int m=0; m<2; m++){
							final int total=p.parseInt(8+4*m)+p.parseInt(9+4*m);
							final int shortCount=id.equals("long256") ? p.parseInt(8+4*m) : Math.min(1, total);
							out.tab().append(p.parseInt(6+4*m)).tab().append(p.parseInt(7+4*m))
								.tab().append(shortCount).tab().append(total-shortCount);
						}
						output.writer.print(out.nl());
					}
				}
			}
		}
	}

	/** ByteBuilder String appends are ASCII-only; the literal alpha panel name
	 * must retain its original UTF-8 bytes in both interchange versions. */
	private static void identity(final ByteBuilder out, final LineParser1 p){
		assert(out!=null && p!=null && p.terms()==14)
			: "Copied window identities come from the complete two-model saved-count schema";
		for(int field=2; field<6; field++){
			p.setBounds(field);
			out.tab().append(p.line(), p.a(), p.b()-p.a());
		}
	}

	/** A selected run must reproduce its own version's full-run rows exactly. */
	private static void compareSelected(final String original, final String changed, final String layout){
		final LineParser1 a=new LineParser1('\t'), b=new LineParser1('\t');
		try(Input left=new Input(original); Input right=new Input(changed)){
			left.reader.nextLine();
			right.reader.nextLine();
			int matched=0;
			for(byte[] row; (row=left.reader.nextLine())!=null;){
				a.set(row);
				if(!a.termEquals(layout, 0)){continue;}
				byte[] other;
				do{
					other=right.reader.nextLine();
					check(other!=null, "Generic fitting must retain every corresponding legacy fit/trial row");
					b.set(other);
				}while(!b.termEquals(layout, 0));
				check(Arrays.equals(row, other), "Selected fitting must preserve every deterministic full-run result byte");
				matched++;
			}
			check(matched>0, "Equivalence comparison must cover real output rows");
			for(byte[] row; (row=right.reader.nextLine())!=null;){
				b.set(row);
				check(!b.termEquals(layout, 0), "Generic fitting cannot append extra selected-layout rows");
			}
		}
	}

	private static void faultCounts(final String source, final String path, final String fault){
		try(Input input=new Input(source); Output output=new Output(path)){
			output.writer.println(input.reader.nextLine());
			final byte[] first=input.reader.nextLine();
			check(first!=null, "Faults start from a complete saved-count fixture");
			if(!fault.equals("missing")){
				final String[] fields=new String(first, StandardCharsets.UTF_8).split("\t");
				if(fault.equals("identity")){fields[3]="999";}
				if(fault.equals("overflow")){fields[5]="65536";}
				if(fault.equals("forbidden")){fields[7]="1";}
				output.writer.println(String.join("\t", fields).getBytes(StandardCharsets.UTF_8));
			}
			for(byte[] row; (row=input.reader.nextLine())!=null;){output.writer.println(row);}
			if(fault.equals("extra")){output.writer.println(first);}
		}
	}

	/** Current rows are synthetic transformations of the same tiny saved fixture;
	 * their purpose is interchange shape and fitting parity, not layout quality. */
	private static void current(final String root, final OptimizeNcrnaPartitionWeights.Config original){
		final OptimizeNcrnaPartitionWeights.Config all=config(original, root+"/v2_fit");
		all.randomDraws=7;
		all.layouts=root+"/v2_layouts.tsv";
		all.counts=root+"/v2_counts.tsv.gz";
		metadataV2(all.layouts);
		countsV2(original.output+".counts.tsv.gz", all.counts);
		final NcrnaDdlLayoutCounts manifest=new NcrnaDdlLayoutCounts(all.layouts, 2);
		check(manifest.schemaVersion==2 && !manifest.historical && manifest.control(all.control)==1,
			"Version2 binds its explicit legacy control independently of row position");
		final int[] terms={2,3,2,1,2};
		final int[] references={0,1,0,0,0};
		final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(all.panels, all.selection, all.models, all.groups);
		try(NcrnaDdlLayoutCounts.Reader reader=new NcrnaDdlLayoutCounts.Reader(all.counts, manifest, data)){
			for(int i=0; i<V2_IDS.length; i++){
				final NcrnaDdlLayoutCounts.Layout layout=manifest.layouts[i];
				check(layout.id.equals(V2_IDS[i]) && layout.terms==terms[i] && layout.referenceTerm==references[i]
					&& !layout.historical && layout.memberOriginPayloadBytes==0
					&& layout.memberKeyPayloadBytes==4L*layout.memberBins,
					"One/two/three-term layouts use only their actual active prefix and plain two-byte model keys");
				check(layout.queryForm.equals(i==4 ? "QUERY_DDL" : "RAW_KEYS"),
					"Both supported query forms retain their explicit metadata identity");
				final char[][] rows=reader.next(true);
				check(rows.length==11,
					"RankSavedNcrnaWindowsTest.Fixture supplies ten alpha windows and one beta window to every layout");
				for(final char[] row:rows){
					check(row.length==8, "Two original models require exactly four interchange slots each");
					for(int model=0; model<2; model++){
						for(int term=terms[i]; term<4; term++){
							check(row[4*model+term]==0, "Inactive model slots must stay zero after reading the current active prefix");
						}
					}
				}
			}
			reader.finish();
		}
		OptimizeNcrnaPartitionWeights.run(all);
		final OptimizeNcrnaPartitionWeights.Config selected=config(all, root+"/v2_cli7");
		check(selected.randomDraws==7, "The fixture configuration copier must preserve a nondefault search budget");
		selected.layouts=all.layouts;
		selected.counts=all.counts;
		selected.fitLayouts=all.control;
		cli(selected, 7);
		for(final String suffix:FIT_SUFFIXES){compareSelected(all.output+suffix, selected.output+suffix, all.control);}
		checkRandomDraws(selected.output+".trials.tsv.gz", all.control, 7);
		checkRandomDraws(all.output+".trials.tsv.gz", "pooled256", 0);
		try(Input input=new Input(selected.output+".comparison.tsv")){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.termEquals(all.control, 0) && p.parseInt(34)==7 && input.reader.nextLine()==null,
				"The actual CLI budget must reach both trial generation and its reported selected-layout summary");
		}
		selected.output=root+"/v2_cli_zero";
		boolean rejected=false;
		try{cli(selected, 0);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(selected.output+".trials.tsv.gz").exists(),
			"randomdraws=0 must fail through the actual CLI before trial publication");
		currentFaults(root, all);
		System.err.println("NCRNA_DDL_V2_ACTIVE_PREFIX_TEST_PASS");
	}

	private static void cli(final OptimizeNcrnaPartitionWeights.Config c, final int draws){
		assert(c!=null && c.fitLayouts!=null) : "The CLI witness must name its selected layout and complete fixture configuration";
		OptimizeNcrnaPartitionWeights.main(new String[]{"layouts="+c.layouts, "counts="+c.counts,
			"control="+c.control, "fitlayouts="+c.fitLayouts, "panels="+c.panels,
			"selection="+c.selection, "ref="+c.references, "identities="+c.identities,
			"ranks="+c.ranks, "out="+c.output, "expectedmodels="+c.models,
			"expectedgroups="+c.groups, "randomdraws="+draws, "t=1"});
	}

	/** Count actual RANDOM rows, rather than trusting a copied budget label. */
	private static void checkRandomDraws(final String path, final String id, final int expected){
		final int[] counts=new int[3], units=new int[3];
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals(id, 0)){continue;}
				final int fit=p.termEquals("A", 2) ? 0 : p.termEquals("B", 2) ? 1 : 2;
				if(p.termEquals("RANDOM", 4)){
					check(p.parseInt(3)>=1 && p.parseInt(3)<=expected, "Random trial ordinals must stay within the actual requested budget");
					counts[fit]++;
				}else if(p.termEquals("UNIT", 4)){
					check(p.parseInt(3)==0, "Every fit starts from the unit-weight control");
					units[fit]++;
				}
			}
		}
		check(units[0]==3 && units[1]==3 && counts[0]==3*expected && counts[1]==3*expected,
			"Each judging direction reports all three populations for exactly the requested random points");
		check(units[2]==0 ? counts[2]==0 : units[2]==3 && counts[2]==3*expected,
			"Any permitted refit must reuse the same initial random population");
	}

	private static void metadataV2(final String path){
		assert(path!=null) : "Synthetic metadata must have a named fixture output, never implicit standard output";
		try(Output out=new Output(path)){
			out.writer.println(NcrnaDdlLayoutCounts.LAYOUT_HEADER_V2);
			for(final String id:V2_IDS){
				final boolean legacy=id.equals(IDS[1]), pooled=id.equals("pooled256");
				final int member=legacy ? 384 : 256, consensus=legacy ? 128 : pooled ? 0 : 256;
				final int slots=member+consensus;
				final ByteBuilder b=new ByteBuilder();
				b.append(2).tab().append(id).append("\tTOY\t7\t10\tMIDDLE_RIGHT\t")
					.append(legacy ? "PARTITIONED" : "MIXED").tab().append(member)
					.tab().append(legacy ? 128 : 0).tab().append(legacy ? 256 : 0)
					.tab().append(pooled ? "POOL_BOTH" : id.equals("long256") ? "LONG_ONLY" : "BOTH")
					.tab().append(consensus).tab().append(pooled ? 5 : 0).append("\t4\t12\tWRAP\t")
					.append(legacy ? "member_short,member_long,consensus_total" : pooled ? "member_total" : "member_total,consensus_total")
					.tab().append(legacy ? 1 : 0).append("\t2").tab().append(4L*member).tab().append(4L*consensus)
					.tab().append(slots).tab().append(slots).append("\t0\t4").tab().append(4L*slots)
					.tab().append(8L*slots).append("\t0\t0\t").append(legacy ? "LEGACY_PARTITIONS" : "LENGTH_TAGGED")
					.tab().append(id.equals("ddl256") ? "QUERY_DDL" : "RAW_KEYS").nl();
				out.writer.print(b);
			}
		}
	}

	/** Preserve the exact control total; other rows are explicitly synthetic
	 * count shapes, not claims about new raw/query-DDL biological performance. */
	private static void countsV2(final String source, final String path){
		assert(source!=null && path!=null && !source.equals(path)) : "Current synthetic counts must preserve their original capture as a separate input";
		final LineParser1 p=new LineParser1('\t');
		try(Output out=new Output(path)){
			out.writer.println(NcrnaDdlLayoutCounts.COUNT_HEADER_V2);
			for(final String id:V2_IDS){
				try(Input input=new Input(source)){
					input.reader.nextLine();
					for(byte[] row; (row=input.reader.nextLine())!=null;){
						p.set(row);
						check(p.terms()==14, "The original two-model capture must retain six identity fields and four counts per model");
						if(!p.termEquals(id.equals("long256") ? "CONSENSUS10_ONLY" : "CONSENSUS_BOTH_LENGTHS", 0)
							|| p.parseInt(1)!=512){continue;}
						final ByteBuilder b=new ByteBuilder().append(id);
						identity(b, p);
						for(int model=0; model<2; model++){
							final int shortCount=p.parseInt(6+4*model), longCount=p.parseInt(7+4*model);
							final int consensus=p.parseInt(8+4*model)+p.parseInt(9+4*model);
							if(id.equals(IDS[1])){b.tab().append(shortCount).tab().append(longCount).tab().append(consensus).tab().append(0);}
							else{b.tab().append(shortCount+longCount).tab().append(id.equals("pooled256") ? 0 : consensus).append("\t0\t0");}
						}
						out.writer.print(b.nl());
					}
				}
			}
		}
	}

	private static void currentFaults(final String root, final OptimizeNcrnaPartitionWeights.Config all){
		final String[] metadataFaults={"shape", "terms", "reference", "query", "control_query", "legacy_geometry", "member_bytes"};
		final int[] columns={-1,16,17,30,30,8,19};
		final String[] values={"extra", "member_total,consensus_short,consensus_long", "1", "UNKNOWN", "QUERY_DDL", "127", "0"};
		final String[] reasons={"versioned field", "score terms", "score terms", "Query form", "no index arrays", "Legacy namespace", "two bytes"};
		for(int i=0; i<metadataFaults.length; i++){
			final String path=root+"/v2_bad_metadata_"+metadataFaults[i]+".tsv";
			mutate(all.layouts, path, 1, i==4 || i==5 ? all.control : V2_IDS[0], columns[i], values[i]);
			boolean failed=false;
			try{new NcrnaDdlLayoutCounts(path, 2).control(all.control);}
			catch(IllegalArgumentException expected){
				check(expected.getMessage().contains(reasons[i]), "Metadata must fail for the intended mismatch: "+expected.getMessage());
				failed=true;
			}
			check(failed, "Invalid current metadata cannot masquerade as a supported plain-key layout");
		}
		final String[] ids={V2_IDS[0], all.control, "pooled256", V2_IDS[0], V2_IDS[0], V2_IDS[0]};
		final int[] slots={7,8,6,5,5,-1};
		final String[] badValues={"1","1","1","-1","65536","extra"};
		for(int i=0; i<ids.length; i++){
			final OptimizeNcrnaPartitionWeights.Config bad=config(all, root+"/v2_rejected_counts_"+i);
			bad.layouts=all.layouts;
			bad.counts=root+"/v2_bad_counts_"+i+".tsv";
			bad.fitLayouts="pooled256";
			mutate(all.counts, bad.counts, 0, ids[i], slots[i], badValues[i]);
			boolean failed=false;
			try{OptimizeNcrnaPartitionWeights.run(bad);}
			catch(IllegalArgumentException expected){
				check(expected.getMessage().contains(i==5 ? "Count rows" : "zero padding"),
					"Count rejection must identify the intended active-prefix/shape fault: "+expected.getMessage());
				failed=true;
			}
			check(failed && !new File(bad.output+".trials.tsv.gz").exists(),
				"Malformed active/padding slots, even in unselected layouts, must fail before fitting");
		}
	}

	/** Alter exactly one named row, preserving all unrelated fixture bytes. */
	private static void mutate(final String source, final String output, final int idColumn,
			final String id, final int column, final String value){
		boolean changed=false;
		try(Input in=new Input(source); Output out=new Output(output)){
			out.writer.println(in.reader.nextLine());
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				final String[] fields=new String(row, StandardCharsets.UTF_8).split("\t", -1);
				if(!changed && fields[idColumn].equals(id)){
					if(column<0){out.writer.println((String.join("\t", fields)+"\t"+value).getBytes(StandardCharsets.UTF_8));}
					else{
						fields[column]=value;
						out.writer.println(String.join("\t", fields).getBytes(StandardCharsets.UTF_8));
					}
					changed=true;
				}else{out.writer.println(row);}
			}
		}
		check(changed, "Every corruption witness must change exactly one existing row: "+id);
	}

	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private static final String[] IDS={"long256", "legacy_member384_consensus128_both", "mixed256"};
	private static final String[] FIT_SUFFIXES={".trials.tsv.gz", ".fits.tsv", ".decision.tsv"};
	private static final String[] V2_IDS={"mixed256", IDS[1], "long256", "pooled256", "ddl256"};
}
