package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Materializes a complete caller population from accepted saved-rank metadata.
 * Derived from template.A_SampleByteFile. No sequences, calling, matching or
 * ranking are performed. Original identities survive new panel-local ordinals.
 * Exact-sequence and same-source labels remain separate, including UNKNOWN.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlCallerPanel{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlCallerPanel.class, false);
		try{
			Shared.setThreads(1);
			final Parser parser=new Parser();
			parser.overwrite=false;
			String mode=null, audits=null;
			int expectedInputs=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("mode")){mode=b;}
				else if(a.equals("audits")){audits=b;}
				else if(a.equals("expectedinputs")){expectedInputs=Parse.parseIntKMG(b);}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown caller-panel parameter: "+arg);}
			}
			require("capture".equals(mode) || "combined".equals(mode), "Declare mode=capture or mode=combined to select an exact saved metadata schema");
			require(parser.in1!=null && parser.out1!=null && audits!=null && !audits.isEmpty()
				&& parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1 && expectedInputs>0 && Shared.threads()==1,
				"Supply complete in=metadata audits=original,additional out=fresh_prefix expectedinputs=N without sampling or overwrite");
			new NcrnaDdlCallerPanel().process(parser.in1, "combined".equals(mode), audits.split(",", -1), parser.out1, expectedInputs);
			pp.outstream.println("NCRNA_DDL_CALLER_PANEL_COMPLETE_METADATA_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final String input, final boolean combined, final String[] audits,
			final String output, final int expectedInputs){
		RrnaResourceIO.fresh(output, ".config", ".inputs.tsv", ".strata.tsv", ".summary.tsv");
		loadInputs(input, combined);
		require(samples.size()==expectedInputs, "Caller populations must retain every declared input, including genomes with no target");
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		for(final String audit:audits){
			require(audit!=null && !audit.isEmpty() && paths.add(audit), "An audit file may be consumed only once");
			loadAudit(audit);
		}
		for(final Sample sample:samples){
			require(sample.labels.size()==sample.targets, "Missing or duplicated exact-sequence labels for original input "+sample.id+": expected="+sample.targets+", actual="+sample.labels.size());
		}
		write(output);
	}

	/** Combined metadata already proves the selected assembly-group partition. */
	private void loadInputs(final String path, final boolean combined){
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final IntHashMap groups=new IntHashMap();
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, combined ? COMBINED_HEADER : CAPTURE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, combined ? 11 : 13);
				final int slot=p.parseInt(combined ? 1 : 0), id=p.parseInt(combined ? 2 : 1);
				final String fasta=text(combined ? 3 : 2), panel=combined ? text(0) : "original";
				final int group=combined ? p.parseInt(4) : -1, targets=p.parseInt(combined ? 5 : 7);
				require(slot>=0 && id>=0 && fasta.startsWith("/") && targets>=0 && !panel.isEmpty()
					&& !ids.contains(id) && paths.add(fasta), "Caller inputs need unique original IDs/absolute paths and nonnegative source slots/target counts");
				if(combined){
					require(group>=0 && !groups.contains(group), "The expanded caller population requires exactly one representative per resolved assembly group");
					long total=0;
					for(int field=6; field<11; field++){
						final int value=p.parseInt(field);
						require(value>=0, "Saved per-input target dispositions cannot be negative");
						total+=value;
					}
					require(total==targets, "Combined metadata must conserve every target, not only rankable loci");
					groups.put(group, samples.size());
				}else{require(slot==samples.size(), "Original capture slots must remain complete and ordered");}
				ids.put(id, samples.size());
				samples.add(new Sample(id, slot, fasta, panel, group, targets));
			}
		}
	}

	/** Consume explicit exact_sequence, never the composite seen_stratum field. */
	private void loadAudit(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, AUDIT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 21);
				final int index=ids.get(p.parseInt(1));
				if(index<0){continue;}
				final Sample sample=samples.get(index);
				require(p.parseInt(0)==sample.slot && sample.fasta.equals(text(2)), "Audit slots/paths must bind to the selected original input, never a new panel ordinal");
				require(p.termEquals("false", 8) && p.termEquals("true", 9), "Saved overlap labels must describe scored target parents, not masks");
				require(sample.group<0 || p.parseInt(15)==sample.group, "Selected assembly groups must agree with each retained audit parent");
				final String parent=text(3);
				require(!parent.isEmpty() && !parent.equals(".") && sample.parents.add(parent), "Each original target parent must receive exactly one label");
				final int exact=state(17), source=state(18);
				final boolean complete=bool(16), common=bool(14);
				require(complete || exact==2, "A missing complete defined sequence must retain UNKNOWN exact membership");
				require(!complete || exact!=2, "Complete defined sequence must have a definite exact-membership result");
				sample.labels.add(new Label(parent, exact, source, complete, common));
			}
		}
	}

	private void write(final String prefix){
		final long[] counts=new long[3];
		try(Output config=new Output(prefix+".config"); Output inputs=new Output(prefix+".inputs.tsv");
			Output strata=new Output(prefix+".strata.tsv")){
			inputs.writer.println(INPUT_HEADER);
			strata.writer.println(STRATA_HEADER);
			for(int slot=0; slot<samples.size(); slot++){
				final Sample sample=samples.get(slot);
				config.writer.print(b.clear().append("in=").append(sample.fastaBytes).nl());
				inputs.writer.print(b.clear().append(slot).tab().append(sample.id).tab().append(sample.fastaBytes).tab()
					.append(sample.panelBytes).tab().append(sample.slot).tab().append(sample.group).tab().append(sample.targets).nl());
				for(final Label label:sample.labels){
					counts[label.exact]++;
					strata.writer.print(b.clear().append(sample.id).tab().append(label.parentBytes).tab().append(STATES[label.exact]).tab()
						.append(STATES[label.source]).tab().append(label.complete).tab().append(label.common).nl());
				}
			}
		}
		try(Output summary=new Output(prefix+".summary.tsv")){
			summary.writer.println("inputs\ttargets\texact_seen\texact_not_seen\texact_unknown");
			summary.writer.print(b.clear().append(samples.size()).tab().append(counts[0]+counts[1]+counts[2])
				.tab().append(counts[0]).tab().append(counts[1]).tab().append(counts[2]).nl());
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Parse Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private void parse(final byte[] row, final int fields){
		p.set(row);
		require(p.terms()==fields, "Saved caller-panel metadata schema changed; expected fields="+fields);
	}
	private String text(final int field){
		p.setBounds(field);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0 && value.indexOf('\r')<0, "Invalid UTF-8 or carriage returns would change a retained caller-input identity/config");
		return value;
	}
	private boolean bool(final int field){
		require(p.termEquals("true", field) || p.termEquals("false", field), "Audit booleans must be explicit true/false");
		return p.termEquals("true", field);
	}
	private int state(final int field){
		for(int i=0; i<STATES.length; i++){if(p.termEquals(STATES[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown exact-sequence/source-locus state");
	}
	private static void require(final boolean valid, final String message){
		assert(valid) : message;
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Sample{
		Sample(final int id_, final int slot_, final String fasta_, final String panel_, final int group_, final int targets_){
			id=id_;
			slot=slot_;
			fasta=fasta_;
			fastaBytes=fasta_.getBytes(StandardCharsets.UTF_8);
			panelBytes=panel_.getBytes(StandardCharsets.UTF_8);
			group=group_;
			targets=targets_;
		}
		final int id, slot, group, targets;
		final String fasta;
		final byte[] fastaBytes, panelBytes;
		final ObjectSet<String> parents=new ObjectSet<String>(String.class);
		final ArrayList<Label> labels=new ArrayList<Label>();
	}
	private static final class Label{
		Label(final String parent_, final int exact_, final int source_, final boolean complete_, final boolean common_){
			parentBytes=parent_.getBytes(StandardCharsets.UTF_8);
			exact=exact_;
			source=source_;
			complete=complete_;
			common=common_;
		}
		final byte[] parentBytes;
		final int exact, source;
		final boolean complete, common;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<Sample> samples=new ArrayList<Sample>();
	private final IntHashMap ids=new IntHashMap();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private static final String[] STATES={"SEEN", "NOT_SEEN", "UNKNOWN"};
	public static final String CAPTURE_HEADER="input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tactual_attempts\ttarget_loci\tseed_cpu_ns\tnative_index_cpu_ns\toracle_cpu_ns\tranking_cpu_ns\tworker_report_cpu_ns";
	public static final String COMBINED_HEADER="panel\tinput_slot\tinput_id\tfasta\tresolved_assembly_group\ttarget_loci\tcommon_loci\tnative_zero_support\tno_compatible_model\tno_window\tambiguous_only";
	public static final String AUDIT_HEADER=NcrnaSavedRankSummary.PARENT_HEADER+"\tcommon_native_rankable\tresolved_assembly_group\tcomplete_defined_query\texact_sequence\tsame_source_locus\tseen_stratum\tmatched_source_locus";
	public static final String INPUT_HEADER="input_slot\tinput_id\tfasta\tsource_panel\tsource_slot\tresolved_assembly_group\ttarget_loci";
	public static final String STRATA_HEADER="input_id\treference_id\texact_sequence\tsame_source_locus\tcomplete_defined_query\tcommon_native_rankable";
}
