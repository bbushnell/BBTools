package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import map.LongIntMap;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;

/** Binds prior decisions to saved input/contig observations without adjudicating them.
 * Versioned readers preserve the original rows, classifications and interval meanings.
 * A current model hit is neither read nor required. All state belongs to one thread;
 * original sequence, annotation and model payloads are never opened.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferencePrior{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	static void process(final String input, final String manifest, final String output,
		final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceProfiles profiles=NcrnaReferenceProfiles.load(input, catalog, roster);
		final NcrnaReferencePrior prior=new NcrnaReferencePrior(catalog, profiles);
		prior.readManifest(manifest, input, output);
		prior.readContigs(input+".contigs.tsv.gz");
		prior.write(output);
	}

	private NcrnaReferencePrior(final NcrnaEvidenceCatalog catalog_, final NcrnaReferenceProfiles profiles_){
		catalog=catalog_;
		profiles=profiles_;
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			final GenomeState state=new GenomeState(genome);
			byId.put(genome.id, state);
			addPath(genome.fasta, state);
			for(int j=0; j<genome.memberships(); j++){addPath(genome.membership(j).suppliedFasta, state);}
		}
	}

	/** Only explicitly recorded paths may identify an input; conflicting claims fail. */
	private void addPath(final String path, final GenomeState state){
		final GenomeState previous=byPath.get(path);
		require(previous==null || previous==state, "Original FASTA path maps to multiple catalog inputs: "+path);
		byPath.put(path, state);
	}

	private void readManifest(final String path, final String input, final String output){
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		final ArrayList<String> paths=new ArrayList<String>();
		paths.add(path);
		paths.add(input+".contigs.tsv.gz");
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Prior manifest needs six fields: "+path);
				final String id=text(p, 0), format=text(p, 1), decisions=text(p, 2), bindings=text(p, 3), authority=text(p, 4);
				final long expected=p.parseLong(5);
				require(!id.isEmpty() && ids.add(id) && !authority.isEmpty() && expected>=0,
					"Prior sources need unique names, explicit authority and nonnegative row counts");
				require(format.equals(PANEL) || format.equals(FMN), "Unknown prior-decision format: "+format);
				require(new File(decisions).isAbsolute() && new File(bindings).isAbsolute(), "Prior source paths must be absolute");
				sources.add(new PriorSource(id, format, decisions, bindings, authority, expected));
				paths.add(decisions);
				paths.add(bindings);
			}
		}
		require(!sources.isEmpty(), "Declare at least one prior source, including its zero-row count when empty");
		require(Tools.testInputFiles(false, true, paths.toArray(new String[paths.size()])), "A declared prior input cannot be read");
		for(final String suffix:OUTPUT_SUFFIXES){paths.add(output+suffix);}
		require(Tools.testForDuplicateFiles(true, paths.toArray(new String[paths.size()])), "Prior source and output paths must be distinct");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Saved contigs        ----------------*/
	/*--------------------------------------------------------------*/

	private void readContigs(final String path){
		final LineParser1 p=new LineParser1('\t');
		GenomeState state=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "input_id\tfasta\tliteral_target\tobserved_length", path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==4, "Saved contig table needs four fields");
				if(state==null || !p.termEquals(state.genome.id, 0)){state=byId.get(text(p, 0));}
				require(state!=null && profiles.selected(state.genome.id), "Saved contig names an input outside the saved population");
				require(same(p, 1, state.encodedPath), "Saved contig FASTA differs from its census input");
				final String target=text(p, 2);
				final int length=p.parseInt(3);
				require(!target.isEmpty() && length>0 && !state.lengths.contains(target), "Literal contigs need unique names and positive lengths within each input");
				state.lengths.put(target, length);
				state.bases=Math.addExact(state.bases, length);
			}
		}
		for(int i=0; i<catalog.size(); i++){
			final GenomeState item=byId.get(catalog.get(i).id);
			if(!profiles.selected(item.genome.id)){continue;}
			require(item.lengths.size()==item.genome.records && item.bases==item.genome.bases,
				"Saved contigs must conserve the selected input's census dimensions: "+item.genome.id);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Original decisions     ----------------*/
	/*--------------------------------------------------------------*/

	private void readPanel(final PriorSource source, final Output out){
		final ObjectMap<String,Binding> bindings=new ObjectMap<String,Binding>(String.class, Binding.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(source.bindings)){
			RrnaResourceIO.header(in.reader, NcrnaCalibrationPrep.ROSTER_HEADER, source.bindings);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11, "Prior panel roster needs eleven fields");
				final String genome=text(p, 1), fasta=text(p, 7);
				require(!genome.isEmpty() && !bindings.contains(genome), "Prior roster genome key must be unique within its source");
				bindings.put(genome, bind(fasta));
			}
		}
		final ObjectMap<String,LocusKeys> keys=new ObjectMap<String,LocusKeys>(String.class, LocusKeys.class);
		try(Input in=new Input(source.decisions)){
			RrnaResourceIO.header(in.reader, CobalaminPanelPrep.HEADER, source.decisions);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				require(p.terms()==10, "Prior panel locus table needs ten fields");
				final String locus=text(p, 0), genome=text(p, 1), target=text(p, 2), family=text(p, 3), assignment=text(p, 4);
				final Binding binding=bindings.get(genome);
				require(binding!=null && !locus.isEmpty(), "Every prior locus must bind its original roster genome");
				boolean known=false;
				for(final String name:CobalaminPanelPrep.FAMILIES){known|=name.equals(family);}
				require(known, "Unknown prior cobalamin source family: "+family);
				final boolean core=assignment.equals(family) || assignment.equals(CobalaminPanelPrep.CLAN_ONLY);
				require(core || assignment.equals(CobalaminPanelPrep.UNRESOLVED) || assignment.equals(CobalaminPanelPrep.UNREVIEWED),
					"Unknown prior cobalamin assignment: "+assignment);
				final long first=p.parseLong(5), last=p.parseLong(6), start=p.parseLong(7), stop=p.parseLong(8);
				require(core ? start>0 && stop>=start : start==0 && stop==0, "Prior assignment and original core-presence fields disagree");
				LocusKeys seen=keys.get(genome);
				if(seen==null){
					seen=new LocusKeys();
					keys.put(genome, seen);
				}
				require(seen.keys.add(family+'\t'+locus), "Repeated original family/locus key within a prior genome: "+genome+" "+locus);
				emit(source, line, binding, locus, genome, target, family, assignment, first, last, start, stop,
					strand(p, 9), "SOURCE_AND_REVIEWED_CORE", "NA", row, out);
			}
		}
	}

	/** Physical-line bindings also match original row fields, preventing stale row-number joins. */
	private void readFmn(final PriorSource source, final Output out){
		final ArrayList<FmnBinding> bindings=new ArrayList<FmnBinding>();
		final LongIntMap byLine=new LongIntMap(16);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(source.bindings)){
			RrnaResourceIO.header(in.reader, FMN_BINDING_HEADER, source.bindings);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "FMN bindings need seven fields");
				final long line=p.parseLong(0);
				require(line>=2 && !byLine.contains(line), "Each FMN decision line needs exactly one original FASTA binding");
				byLine.put(line, bindings.size());
				bindings.add(new FmnBinding(bind(text(p, 1)), text(p, 2), p.parseLong(3), p.parseLong(4), strand(p, 5), text(p, 6)));
			}
		}
		try(Input in=new Input(source.decisions)){
			RrnaResourceIO.header(in.reader, FMN_HEADER, source.decisions);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				require(p.terms()==16 && byLine.contains(line), "Every FMN decision needs its exact source schema and FASTA binding");
				final FmnBinding binding=bindings.get(byLine.get(line));
				final String target=text(p, 0), assignment=text(p, 4);
				final long first=p.parseLong(1), last=p.parseLong(2);
				final char strand=strand(p, 3);
				require(!binding.used && binding.target.equals(target) && binding.first==first && binding.last==last
					&& binding.strand==strand && binding.assignment.equals(assignment), "FMN binding fields differ from their original decision row");
				require(assignment.equals("ACCEPTED_NOVEL_FMN") || assignment.equals("UNRESOLVED"), "Unknown prior FMN classification");
				binding.used=true;
				emit(source, line, binding.binding, "NA", "NA", target, "RF00050", assignment, first, last, 0, 0,
					strand, "PRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS", text(p, 12), row, out);
			}
		}
		for(final FmnBinding binding:bindings){require(binding.used, "FMN binding has no original decision row");}
	}

	private Binding bind(final String fasta){
		final GenomeState state=byPath.get(fasta);
		require(state!=null && profiles.selected(state.genome.id), "Original FASTA is absent from the saved population: "+fasta);
		return new Binding(fasta, state);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output             ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		try(Output rows=new Output(prefix+OUTPUT_SUFFIXES[0])){
			rows.writer.println(ROW_HEADER);
			for(final PriorSource source:sources){
				if(source.format.equals(PANEL)){readPanel(source, rows);}
				else{readFmn(source, rows);}
				require(source.rows==source.expected, "Prior row count changed for "+source.id+": "+source.rows+" expected "+source.expected);
			}
		}
		try(Output out=new Output(prefix+OUTPUT_SUFFIXES[1])){
			out.writer.println("source_id\tassignment\tdecision_rows");
			for(final PriorSource source:sources){
				long sum=0;
				for(final Count count:source.orderedCounts){
					field(buffer.clear(), source.id).tab();
					field(buffer, count.assignment).tab().append(count.rows).nl();
					out.writer.print(buffer);
					sum=Math.addExact(sum, count.rows);
				}
				require(sum==source.rows, "Every original decision contributes once to its unchanged assignment");
			}
		}
		// Completion follows every reader, binding check and output closure.
		try(Output out=new Output(prefix+OUTPUT_SUFFIXES[2])){
			out.writer.println("sources\tdecision_rows\tinputs_observed\tcatalog_inputs\tevidence_status\tstatus");
			buffer.clear().append(sources.size()).tab().append(rowsWritten).tab().append(profiles.observed()).tab().append(catalog.size()).tab();
			field(buffer, profiles.evidenceState()).append("\tPRIOR_DECISIONS_BOUND_NOT_NEWLY_ADJUDICATED\n");
			out.writer.print(buffer);
		}
	}

	private void emit(final PriorSource source, final long line, final Binding binding, final String locus,
		final String genome, final String target, final String family, final String assignment,
		final long first, final long last, final long start, final long stop, final char strand,
		final String semantics, final String rowAuthority, final byte[] original, final Output out){
		assert(binding!=null && original!=null && line>=2) : "Every output must retain a physical source row and its original FASTA binding";
		field(buffer.clear(), target);
		// ByteBuilder.toString is ASCII-only; saved contig keys retain UTF-8.
		final String encodedTarget=new String(buffer.array, 0, buffer.length(), StandardCharsets.UTF_8);
		final int length=binding.state.lengths.get(encodedTarget);
		require(length>0, "Prior decision has no observed literal contig in its input: "+binding.fasta+" "+target);
		require(first>0 && last>=first && last<=length && (start==0 && stop==0 || start>0 && stop>=start && stop<=length),
			"Prior source/core interval leaves its original contig: "+source.id+" line "+line);
		field(buffer.clear(), source.id).tab();
		field(buffer, source.decisions).tab().append(line).tab();
		field(buffer, locus).tab();
		field(buffer, genome).tab();
		field(buffer, binding.fasta).tab();
		field(buffer, binding.state.genome.id).tab();
		field(buffer, binding.state.genome.fasta).tab();
		field(buffer, target).tab().append(length).tab();
		field(buffer, family).tab();
		field(buffer, assignment).tab().append(first).tab().append(last).tab().append(start).tab().append(stop).tab().append(strand).tab();
		field(buffer, semantics).tab();
		field(buffer, source.authority).tab();
		field(buffer, rowAuthority).tab();
		BuildNcrnaReference.field(buffer, original, 0, original.length).nl();
		out.writer.print(buffer);
		Count count=source.counts.get(assignment);
		if(count==null){
			count=new Count(assignment);
			source.counts.put(assignment, count);
			source.orderedCounts.add(count);
		}
		count.rows++;
		source.rows++;
		rowsWritten++;
	}

	/*--------------------------------------------------------------*/
	/*----------------             State             ----------------*/
	/*--------------------------------------------------------------*/

	private static final class GenomeState{
		GenomeState(final Genome genome_){
			genome=genome_;
			encodedPath=field(new ByteBuilder(), genome.fasta).toBytes();
		}
		final Genome genome;
		final byte[] encodedPath;
		final ObjectIntMap<String> lengths=new ObjectIntMap<String>(String.class);
		long bases;
	}
	private static final class Binding{
		Binding(final String fasta_, final GenomeState state_){
			fasta=fasta_;
			state=state_;
		}
		final String fasta;
		final GenomeState state;
	}
	private static final class LocusKeys{
		final ObjectSet<String> keys=new ObjectSet<String>(String.class);
	}
	private static final class FmnBinding{
		FmnBinding(final Binding binding_, final String target_, final long first_, final long last_, final char strand_, final String assignment_){
			binding=binding_;
			target=target_;
			first=first_;
			last=last_;
			strand=strand_;
			assignment=assignment_;
		}
		final Binding binding;
		final String target, assignment;
		final long first, last;
		final char strand;
		boolean used;
	}
	private static final class PriorSource{
		PriorSource(final String id_, final String format_, final String decisions_, final String bindings_, final String authority_, final long expected_){
			id=id_;
			format=format_;
			decisions=decisions_;
			bindings=bindings_;
			authority=authority_;
			expected=expected_;
		}
		final String id, format, decisions, bindings, authority;
		final long expected;
		long rows;
		final ObjectMap<String,Count> counts=new ObjectMap<String,Count>(String.class, Count.class);
		final ArrayList<Count> orderedCounts=new ArrayList<Count>();
	}
	private static final class Count{
		Count(final String assignment_){assignment=assignment_;}
		final String assignment;
		long rows;
	}
	private static char strand(final LineParser1 p, final int term){
		require(p.termEquals("+", term) || p.termEquals("-", term), "Prior strand must remain literal + or -");
		return p.termEquals("+", term) ? '+' : '-';
	}
	private static String text(final LineParser1 p, final int term){
		p.setBounds(term);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static boolean same(final LineParser1 p, final int term, final byte[] bytes){
		p.setBounds(term);
		if(p.b()-p.a()!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(p.line()[p.a()+i]!=bytes[i]){return false;}}
		return true;
	}
	private static ByteBuilder field(final ByteBuilder b, final String value){return BuildNcrnaReference.field(b, value);}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private final NcrnaEvidenceCatalog catalog;
	private final NcrnaReferenceProfiles profiles;
	private final ObjectMap<String,GenomeState> byId=new ObjectMap<String,GenomeState>(String.class, GenomeState.class);
	private final ObjectMap<String,GenomeState> byPath=new ObjectMap<String,GenomeState>(String.class, GenomeState.class);
	private final ArrayList<PriorSource> sources=new ArrayList<PriorSource>();
	private final ByteBuilder buffer=new ByteBuilder();
	private long rowsWritten;
	static final String[] OUTPUT_SUFFIXES={".prior_decisions.tsv.gz", ".prior_counts.tsv", ".prior_summary.tsv"};
	static final String MANIFEST_HEADER="source_id\tformat\tdecisions\tbindings\tauthority\texpected_rows";
	static final String FMN_BINDING_HEADER="decision_line\tfasta\tseqid\tstart1\tstop1\tstrand\tclassification";
	static final String FMN_HEADER="seqid\tstart1\tstop1\tstrand\tclassification\twinning_model\tverifier_identity\tmodel_length\toriginal_annotated_members\tbest_CM_bits\tCM_GA_bits\tCM_status\tdecision_source\tindependent_reviewer\ttrace_packet_sha80\tCM_packet_sha80";
	static final String ROW_HEADER="source_id\tdecision_file\tdecision_line\toriginal_locus\toriginal_genome\toriginal_fasta\tinput_id\tresolved_fasta\tliteral_target\tcontig_length\tsource_family\tassignment\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand\tinterval_semantics\tsource_authority\trow_authority\toriginal_row";
	private static final String PANEL="cobalamin_panel_v1", FMN="fmn_adjudication_v1";
}
