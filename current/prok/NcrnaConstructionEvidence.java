package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;

import map.IntHashMap;
import map.LongHashSet;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/** Joins original FMN members to saved source/native-hit correspondence.
 * Cohorts retain separate input/event namespaces. Original source-line witnesses
 * select links; direct bound-input/contig/core overlaps also recover evidence
 * where the saved catalog had no source GFF. Link bases remain explicit, and
 * measured core overlap, not a broad source span, quantifies them.
 * Other-family hits remain evidence for review, never automatic deletion. All
 * RF00050 events are retained separately, including partial and rejected events.
 * This performs no search, sequence extraction, alias proof or admission.
 * @author Ganyu
 */
public final class NcrnaConstructionEvidence{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionEvidence.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String bindings=null, evaluation=null, discovery=null;
			int expected=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("bindings")){bindings=b;}
				else if(a.equals("evaluation")){evaluation=b;}
				else if(a.equals("discovery")){discovery=b;}
				else if(a.equals("expectedmembers")){expected=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown construction evidence parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && bindings!=null && evaluation!=null && discovery!=null && parser.out1!=null && expected>0,
				"Use -ea and supply bindings=prefix evaluation=tier_directory discovery=tier_directory expectedmembers=N out=fresh_prefix");
			require(parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Saved evidence must not be sampled or overwritten");
			process(bindings, evaluation, discovery, parser.out1, expected);
			pp.outstream.println("NCRNA_CONSTRUCTION_EVIDENCE_COMPLETE_NOT_ADMISSION");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void process(final String bindings, final String evaluation, final String discovery,
			final String output, final int expected){
		require(bindings!=null && evaluation!=null && discovery!=null && output!=null && expected>0,
			"Complete binding and both evidence cohorts are required");
		RrnaResourceIO.fresh(output, ".old_evidence.tsv.gz", ".links.tsv.gz", ".fmn_events.tsv.gz", ".summary.tsv");
		final NcrnaConstructionEvidence report=new NcrnaConstructionEvidence();
		report.members(bindings+".members.tsv.gz", expected);
		report.bindings(bindings+".bindings.tsv");
		try(Output old=new Output(output+".old_evidence.tsv.gz"); Output links=new Output(output+".links.tsv.gz");
				Output fmn=new Output(output+".fmn_events.tsv.gz"); Output summary=new Output(output+".summary.tsv")){
			old.writer.println("cohort\told_member_id\tinput_id\tsource_rows_expected\tsource_rows_observed\tdistinct_linked_events\tfull_GA_FMN\tfull_GA_other_family\tmax_other_core_overlap_bp\told_core_length\tclassification\teligibility_status");
			links.writer.println("cohort\told_member_id\told_core_length\tcore_overlap_bp\tlink_basis\t"+EVENT_HEADER);
			fmn.writer.println("cohort\tnative_included_full_model\t"+EVENT_HEADER);
			summary.writer.println("cohort\told_member_bindings\tmissing_source_bindings\tdistinct_member_event_links\tall_native_events\tFMN_native_events\tnative_full_GA_FMN_events\tstatus");
			report.cohort(report.evaluation, evaluation, old, links, fmn, summary);
			report.cohort(report.discovery, discovery, old, links, fmn, summary);
		}
	}

	private NcrnaConstructionEvidence(){}

	/*--------------------------------------------------------------*/
	/*----------------       Original Identity      ----------------*/
	/*--------------------------------------------------------------*/

	private void members(final String path, final int expected){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, MEMBER_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 32);
				require(p.termEquals("true", 0) && p.termEquals("UNVERIFIED", 3) && p.termEquals("NOT_EVALUATED", 4),
					"This join consumes original identity bindings, not a silently substituted eligibility ledger");
				final String id=p.parseString(5), genome=p.parseString(6);
				require(id.startsWith("RF00050_") && ids.add(id), "Original member IDs must remain unique FMN identities");
				final String contig=p.parseString(11);
				final char strand=strand(12);
				final int length=p.parseInt(14), limit=p.parseInt(20);
				p.setBounds(13);
				final int[] intervals=intervals(row, p.a(), p.b());
				long total=0;
				int previous=0;
				for(int i=0; i<intervals.length; i+=2){
					final int start=intervals[i], stop=intervals[i+1];
					require(start>previous && stop>=start && stop<=limit, "Bound core intervals must be disjoint and inside their contig");
					total+=stop-start+1L;
					previous=stop;
				}
				require(length>0 && length==total && !contig.isEmpty(), "Bound core length must conserve its actual intervals");
				p.setBounds(25);
				final long[] lines=sourceLines(row, p.a(), p.b());
				final LongHashSet seen=new LongHashSet();
				for(int i=0; i<lines.length; i++){
					require(lines[i]>0 && seen.add(lines[i]), "Original source-line witnesses must be positive and unique");
				}
				final Member member=new Member(id, contig, strand, length, intervals, lines);
				MemberGroup group=genomes.get(genome);
				if(group==null){
					group=new MemberGroup();
					genomes.put(genome, group);
				}
				group.members.add(member);
				members.add(member);
			}
		}
		require(members.size()==expected, "Original member binding count must match its declared complete population");
	}

	private void bindings(final String path){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BINDING_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 6);
				final String genome=p.parseString(0), input=p.parseString(4), gff=p.parseString(2);
				require(p.parseInt(4)>0, "Original input IDs must be positive");
				final Cohort cohort=p.termEquals("evaluation", 3) ? evaluation : p.termEquals("discovery", 3) ? discovery : null;
				require(cohort!=null && gff.startsWith("/") && p.termEquals("LITERAL_CATALOG_PATH_NOT_SEQUENCE_PROOF", 5)
					&& seen.add(genome+'\t'+cohort.name+'\t'+input), "Original cohort/input bindings must remain unique and explicitly path-based");
				final MemberGroup group=genomes.get(genome);
				if(group==null){continue;}
				final byte[] gffBytes=gff.getBytes(java.nio.charset.StandardCharsets.UTF_8);
				for(final Member member:group.members){
					final State state=new State(member, input);
					cohort.states.add(state);
					member.bindings++;
					for(final long line:member.lines){
						key.clear().append(input).tab().append(gffBytes).tab().append(line);
						probe.set(key.array, 0, key.length());
						Source source=cohort.sources.get(probe);
						if(source==null){
							source=new Source();
							cohort.sources.put(probe.copy(), source);
						}
						source.states.add(state);
					}
				}
			}
		}
		for(final Member member:members){require(member.bindings>0, "Original member has no cohort/input binding: "+member.id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Evidence         ----------------*/
	/*--------------------------------------------------------------*/

	private void cohort(final Cohort cohort, final String directory, final Output old,
			final Output linksOut, final Output fmnOut, final Output summary){
		readSources(cohort, new File(directory, "source_provenance.tsv.gz").getPath());
		indexUnwitnessed(cohort);
		readLinks(cohort, new File(directory, "source_event_links.tsv.gz").getPath());
		long expectedEvents=0, expectedFmn=-1;
		final ObjectSet<String> families=new ObjectSet<String>(String.class);
		try(Input in=new Input(new File(directory, "family_counts.tsv").getPath())){
			RrnaResourceIO.header(in.reader, FAMILY_HEADER, directory);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 12);
				final long count=p.parseLong(1);
				require(count>=0 && families.add(p.parseString(0)), "Saved raw-event family counts must be nonnegative and unique");
				expectedEvents=Math.addExact(expectedEvents, count);
				if(p.termEquals("FMN", 0)){expectedFmn=count;}
			}
		}
		require(expectedFmn>=0, "The complete cohort must explicitly report FMN, including zero events");
		long rows=0, fmn=0, fullFmn=0;
		try(Input in=new Input(new File(directory, "event_provenance.tsv.gz").getPath())){
			RrnaResourceIO.header(in.reader, EVENT_HEADER, directory);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 30);
				final int event=p.parseInt(0);
				require(event>0 && event==++rows, "Native event IDs must retain the complete consecutive cohort sequence");
				final boolean own=p.termEquals("RF00050", 6);
				require(own==p.termEquals("FMN", 7), "FMN family and RF00050 model identity must agree");
				// Bind by the actual catalog input and literal contig before applying
				// strand/core geometry. This is independent of captured GFF availability.
				MemberStates site=null;
				if(cohort.coordinateBindings>0){
					p.setBounds(1);
					final int inputStart=p.a();
					p.setBounds(2);
					probe.set(row, inputStart, p.b());
					site=cohort.sites.get(probe);
				}
				if(site!=null && !p.termEquals("INVALID_INTERVAL", 22)){
					final int start=p.parseInt(3), stop=p.parseInt(4);
					final char strand=strand(5);
					require(start>0 && stop>=start, "Coordinate evidence requires a valid native interval");
					for(final State state:site.states){
						if(strand==state.member.strand && state.member.overlap(start, stop)>0){link(cohort, event, state, 2);}
					}
				}
				final int index=cohort.eventIndex.get(event);
				if(!own && index<0){continue;}
				final boolean included=bool(17);
				final int from=p.parseInt(8), to=p.parseInt(9), modelLength=p.parseInt(10), truncation=p.parseInt(11);
				require(from>0 && to>=from && to<=modelLength && truncation>=0 && truncation<=3, "Saved model span/truncation is invalid");
				final boolean full=included && from==1 && to==modelLength && truncation==0;
				if(own){
					fmn++;
					if(full){fullFmn++;}
					fmnOut.writer.print(buffer.clear().append(cohort.name).tab().append(full).tab().append(row).nl());
				}
				if(index<0){continue;}
				final EventLinks linked=cohort.links.get(index);
				linked.observed=true;
				final int start=p.parseInt(3), stop=p.parseInt(4);
				final char strand=strand(5);
				require(start>0 && stop>=start && !p.termEquals("INVALID_INTERVAL", 22), "A valid source link cannot point to an invalid native interval");
				for(int j=0; j<linked.states.size(); j++){
					final State state=linked.states.get(j);
					final Member member=state.member;
					require(p.termEquals(state.input, 1) && p.termEquals(member.contig, 2) && strand==member.strand,
						"Native links must retain their input, literal contig and strand, not merely matching event numbers");
					final long overlap=member.overlap(start, stop);
					state.events++;
					if(overlap>0 && full){
						if(own){state.ownFull++;}
						else{
							state.otherFull++;
							state.maxOther=Math.max(state.maxOther, overlap);
						}
					}
					linksOut.writer.print(buffer.clear().append(cohort.name).tab().append(member.id).tab().append(member.length)
						.tab().append(overlap).tab().append(LINK_BASES[linked.bases.get(j)]).tab().append(row).nl());
				}
			}
		}
		require(rows==expectedEvents && fmn==expectedFmn, "Complete native events must conserve the saved whole-cohort and FMN counts");
		for(final EventLinks linked:cohort.links){require(linked.observed, "A saved member link names an absent native event");}
		long missing=0, associations=0;
		for(final State state:cohort.states){
			final boolean absent=state.observedSources!=state.member.lines.length;
			if(absent){missing++;}
			associations+=state.events;
			final String status=absent ? "MISSING_SOURCE_WITNESS" : state.ownFull>0 && state.otherFull>0 ? "MIXED_FULL_MODEL_EVIDENCE"
				: state.ownFull>0 ? "OWN_FULL_FMN" : state.otherFull>0 ? "OTHER_FULL_MODEL_EVIDENCE" : "NO_FULL_GA_SUPPORT";
			old.writer.print(buffer.clear().append(cohort.name).tab().append(state.member.id).tab().append(state.input).tab()
				.append(state.member.lines.length).tab().append(state.observedSources).tab().append(state.events).tab()
				.append(state.ownFull).tab().append(state.otherFull).tab().append(state.maxOther).tab().append(state.member.length)
				.tab().append(status).append("\tNOT_EVALUATED\n"));
		}
		summary.writer.print(buffer.clear().append(cohort.name).tab().append(cohort.states.size()).tab().append(missing).tab()
			.append(associations).tab().append(rows).tab().append(fmn).tab().append(fullFmn).append("\tEVIDENCE_JOIN_NOT_ADMISSION_OR_ALIAS_PROOF\n"));
	}

	private void readSources(final Cohort cohort, final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, SOURCE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 13);
				final Source source=source(cohort, row);
				if(source==null){continue;}
				require(!source.observed, "An original input/GFF/line source witness is duplicated");
				source.observed=true;
				for(final State state:source.states){
					require(p.termEquals(state.member.contig, 3) && strand(7)==state.member.strand
						&& state.member.overlap(p.parseInt(5), p.parseInt(6))>0, "Saved source witness no longer overlaps its original member on the original strand");
					state.observedSources++;
				}
			}
		}
	}

	/** Build the fallback index only for wholly or partly missing source witnesses. */
	private void indexUnwitnessed(final Cohort cohort){
		assert(cohort.coordinateBindings==0) : "Source witnesses must be counted before indexing unresolved bindings once";
		for(final State state:cohort.states){
			if(state.observedSources==state.member.lines.length){continue;}
			key.clear().append(state.input).tab().append(state.member.contig.getBytes(java.nio.charset.StandardCharsets.UTF_8));
			probe.set(key.array, 0, key.length());
			MemberStates site=cohort.sites.get(probe);
			if(site==null){
				site=new MemberStates();
				cohort.sites.put(probe.copy(), site);
			}
			site.states.add(state);
			cohort.coordinateBindings++;
		}
	}

	private void readLinks(final Cohort cohort, final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, LINK_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 7);
				final Source source=source(cohort, row);
				if(source==null){continue;}
				require(source.observed, "An event link cannot replace a missing original source witness");
				if(!bool(4)){continue;}
				final int event=p.parseInt(3);
				require(event>0, "Native linked event IDs must be positive");
				for(final State state:source.states){link(cohort, event, state, 1);}
			}
		}
	}

	/** One retained member/event association, with all independently observed bases. */
	private static void link(final Cohort cohort, final int event, final State state, final int basis){
		assert(event>0 && (basis==1 || basis==2)) : "Saved source and bound-coordinate links share one positive native event namespace";
		int index=cohort.eventIndex.get(event);
		if(index<0){
			index=cohort.links.size();
			cohort.eventIndex.put(event, index);
			cohort.links.add(new EventLinks());
		}
		final EventLinks links=cohort.links.get(index);
		for(int i=0; i<links.states.size(); i++){
			if(links.states.get(i)==state){
				links.bases.set(i, links.bases.get(i)|basis);
				return;
			}
		}
		links.states.add(state);
		links.bases.add(basis);
	}

	private Source source(final Cohort cohort, final byte[] row){
		p.setBounds(2);
		probe.set(row, 0, p.b());
		return cohort.sources.get(probe);
	}
	private void parse(final byte[] row, final int terms){
		p.set(row);
		if(p.terms()!=terms){throw new IllegalArgumentException("Unexpected saved construction/evidence schema: expected fields="+terms);}
	}
	private char strand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "A bound construction strand must be explicit");
		return p.termEquals("+", column) ? '+' : '-';
	}
	private boolean bool(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Saved booleans must be literal true/false");
		return p.termEquals("true", column);
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/** Retain primitive coordinates without allocating transient strings per interval. */
	private static int[] intervals(final byte[] row, final int begin, final int end){
		require(begin<end, "Bound core intervals cannot be empty");
		final IntList list=new IntList();
		int start=begin;
		for(int i=begin; i<=end; i++){
			if(i==end || row[i]==','){
				int dash=start;
				while(dash<i && row[dash]!='-'){dash++;}
				require(dash>start && dash<i-1, "Bound core intervals require start-stop pairs");
				final long left=positive(row, start, dash), right=positive(row, dash+1, i);
				require(left<=Integer.MAX_VALUE && right<=Integer.MAX_VALUE, "Core coordinates exceed supported contig length");
				list.add((int)left);
				list.add((int)right);
				start=i+1;
			}
		}
		return list.toArray();
	}
	private static long[] sourceLines(final byte[] row, final int begin, final int end){
		require(begin<end, "Original source witnesses cannot be empty");
		final LongList list=new LongList();
		int start=begin;
		for(int i=begin; i<=end; i++){
			if(i==end || row[i]==','){
				list.add(positive(row, start, i));
				start=i+1;
			}
		}
		return list.toArray();
	}
	private static long positive(final byte[] row, final int begin, final int end){
		require(begin<end, "Original coordinates and source lines must contain decimal digits");
		long value=0;
		for(int i=begin; i<end; i++){
			final int digit=row[i]-'0';
			require(digit>=0 && digit<=9 && value<=(Long.MAX_VALUE-digit)/10,
				"Original coordinates and source lines require nonoverflowing positive decimals");
			value=10*value+digit;
		}
		require(value>0, "Original coordinates and source lines must be positive");
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Member{
		Member(final String id_, final String contig_, final char strand_, final int length_, final int[] intervals_, final long[] lines_){
			id=id_;
			contig=contig_;
			strand=strand_;
			length=length_;
			intervals=intervals_;
			lines=lines_;
		}
		long overlap(final int start, final int stop){
			long result=0;
			for(int i=0; i<intervals.length; i+=2){result+=Math.max(0L, (long)Math.min(stop, intervals[i+1])-Math.max(start, intervals[i])+1);}
			assert(result<=length) : "Core overlap sums disjoint retained intervals, never the outer span";
			return result;
		}
		final String id, contig;
		final char strand;
		final int length;
		final int[] intervals;
		final long[] lines;
		int bindings;
	}
	private static final class MemberGroup{final ArrayList<Member> members=new ArrayList<Member>();}
	private static final class State{
		State(final Member member_, final String input_){
			member=member_;
			input=input_;
		}
		final Member member;
		final String input;
		int observedSources, events, ownFull, otherFull;
		long maxOther;
	}
	private static final class Source{
		final ArrayList<State> states=new ArrayList<State>();
		boolean observed;
	}
	private static final class EventLinks{
		final ArrayList<State> states=new ArrayList<State>();
		final IntList bases=new IntList();
		boolean observed;
	}
	private static final class MemberStates{final ArrayList<State> states=new ArrayList<State>();}
	private static final class Cohort{
		Cohort(final String name_){name=name_;}
		final String name;
		final ObjectMap<Target, Source> sources=new ObjectMap<Target, Source>(Target.class, Source.class);
		final ObjectMap<Target, MemberStates> sites=new ObjectMap<Target, MemberStates>(Target.class, MemberStates.class);
		final ArrayList<State> states=new ArrayList<State>();
		final ArrayList<EventLinks> links=new ArrayList<EventLinks>();
		final IntHashMap eventIndex=new IntHashMap();
		int coordinateBindings;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder key=new ByteBuilder(), buffer=new ByteBuilder();
	private final Target probe=new Target();
	private final ObjectMap<String, MemberGroup> genomes=new ObjectMap<String, MemberGroup>(String.class, MemberGroup.class);
	private final ArrayList<Member> members=new ArrayList<Member>();
	private final Cohort evaluation=new Cohort("evaluation"), discovery=new Cohort("discovery");

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] LINK_BASES={"", "SOURCE_LINE", "BOUND_CORE_COORDINATE", "SOURCE_LINE_AND_BOUND_CORE"};
	public static final String MEMBER_HEADER="old_seed_input\told_consensus_member\tidentity_holdout\tsequence_alias_status\teligibility_status\t"+NcrnaConstructionBindings.LOCI_HEADER
		+"\tassignment_full_record_id\tassignment_source_tid\tassignment_group\tassignment_cluster\tassignment_model\tassignment_status";
	public static final String BINDING_HEADER="original_genome\toriginal_fasta\toriginal_gff\tcohort\tinput_id\tbasis";
	public static final String SOURCE_HEADER="input_id\tgff\tsource_line\ttarget\tsource_type\tstart\tstop\tstrand\tdeclared_family\tselection_bits\tdisposition\tsource_flags\tqualification";
	public static final String LINK_HEADER="input_id\tgff\tsource_line\tevent_id\tsame_strand_valid\tfamily_agreement\tevent_tier";
	public static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
	public static final String FAMILY_HEADER="family\traw_events\tCONFIDENT\tMARGINAL\tBELOW_MARGINAL\tUNCERTAIN_FLOOR\tREJECTED_WRONG_FAMILY\tCOVERED_BY_PROTECTED_PRIOR\tINVALID_INTERVAL\tscored_annotations\tmasked_annotations\tprotected_priors";
}
