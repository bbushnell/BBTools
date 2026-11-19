package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.function.Consumer;

import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.ListNum;

/** Exact-sequence and corrected source-locus overlap of frozen kept voters.
 * Reads saved windows and construction inputs only. No caller or alignment.
 * Source-locus presence requires a resolved assembly, identical contig and strand,
 * and strictly greater than half overlap of both corrected source and query spans.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlLocusOverlap{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlLocusOverlap.class, false);
		try{
			final Config config=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("members")){config.members=b;}
				else if(a.equals("ledger")){config.ledger=b;}
				else if(a.equals("assignments")){config.assignments=b;}
				else if(a.equals("census")){config.census=b;}
				else if(a.equals("summary")){config.summary=b;}
				else if(a.equals("ref")){config.references=b;}
				else if(a.equals("genomes")){config.genomes=b;}
				else if(a.equals("groups")){config.groups=b;}
				else if(a.equals("assemblies")){config.assemblies=b;}
				else if(a.equals("proxy")){config.proxy=b;}
				else if(a.equals("recovered")){config.recovered=b;}
				else if(a.equals("overlays")){config.overlays=b.split(",", -1);}
				else if(a.equals("expectedinputs")){config.inputs=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgenomes")){config.genomeCount=Parse.parseIntKMG(b);}
				else if(a.equals("expectedtargets")){config.targets=Parse.parseIntKMG(b);}
				else if(a.equals("expectedcommon")){config.common=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmodels")){config.models=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmethods")){config.methods=Parse.parseIntKMG(b);}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown locus-overlap parameter: "+arg);}
			}
			require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use in=saved_rank_prefix out=fresh_prefix without overwrite, append or sampling");
			config.input=parser.in1;
			config.output=parser.out1;
			process(config);
			pp.outstream.println("NCRNA_DDL_LOCUS_OVERLAP_COMPLETE_SAVED_WINDOWS_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void process(final Config config){
		require(config!=null && config.inputs>0 && config.targets>0 && config.common>=0
			&& config.models>0 && config.methods>0 && config.genomeCount>0 && config.overlays!=null
			&& config.input!=null && config.output!=null && config.members!=null && config.ledger!=null
			&& config.assignments!=null && config.census!=null && config.summary!=null && config.references!=null
			&& config.genomes!=null && config.groups!=null && config.assemblies!=null && config.proxy!=null && config.recovered!=null,
			"Declare complete positive input/model/method/target counts and the corrected source overlays");
		RrnaResourceIO.fresh(config.output, ".loci.tsv", ".labels.tsv", ".summary.tsv");
		Shared.TRIM_READ_DESCRIPTION=false;
		final NcrnaDdlLocusOverlap audit=new NcrnaDdlLocusOverlap(config);
		audit.bindings();
		audit.members();
		audit.positions(config.proxy, 0);
		audit.positions(config.recovered, 1);
		for(final String overlay:config.overlays){audit.positions(overlay, 2);}
		audit.parents();
		audit.common();
		audit.windows();
		audit.write();
	}

	private NcrnaDdlLocusOverlap(final Config config){
		c=config;
		samples=new Sample[c.inputs];
	}

	/*--------------------------------------------------------------*/
	/*----------------       Input Bindings         ----------------*/
	/*--------------------------------------------------------------*/

	private void bindings(){
		final ObjectIntMap<String> paths=new ObjectIntMap<String>(String.class);
		try(Input in=new Input(c.groups)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GROUP_HEADER, c.groups);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 14);
				final int group=p.parseInt(10);
				require(group>=0, "Resolved assembly groups must be nonnegative");
				require(p.termEquals("true", 12) || p.termEquals("false", 12), "Unresolved assembly nominations must remain explicit");
				if(p.termEquals("true", 12)){uncertainGroups.put(group, 1);}
				for(int column=2; column<=3; column++){
					final String path=text(column);
					require(path.startsWith("/") && (!paths.contains(path) || paths.get(path)==group), "A literal or canonical source path cannot change resolved assembly");
					paths.put(path, group);
				}
			}
		}
		try(Input in=new Input(c.genomes)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GENOME_HEADER, c.genomes);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				final String source=text(1), path=text(7);
				require(p.parseInt(0)==genomes.size() && !genomeIds.contains(source), "Source genomes retain unique names and complete ordinal order");
				genomeIds.put(source, genomes.size());
				genomes.add(new Genome(source, path, paths.get(path)));
			}
		}
		require(genomes.size()==c.genomeCount, "Source binding must retain the complete declared original genome roster");
		try(Input in=new Input(c.assemblies)){
			RrnaResourceIO.header(in.reader, ASSEMBLY_HEADER, c.assemblies);
			int slot=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 7);
				require(slot<samples.length && p.parseInt(0)==slot, "Every sample must retain its original assembly-census slot");
				final String path=text(2);
				final int group=p.parseInt(3);
				require(group>=-1 && (group<0 || paths.get(path)==group), "The sample assembly result must agree with its accepted path binding");
				samples[slot++]=new Sample(p.parseInt(1), path, group);
			}
			require(slot==samples.length, "Keep every sample input, including those without a target locus");
		}
	}

	private void members(){
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> references=new ArrayList<byte[]>();
		sequences(c.references, read->{
			require(read.bases!=null && read.length()>0, "Every current model needs its complete reference sequence");
			names.add(RrnaResourceIO.first(read.id));
			references.add(read.bases);
		});
		require(names.size()==c.models, "Current reference order must contain the declared complete model set");
		checkedMembers=NcrnaDdlMembers.load(c.members, c.ledger, c.assignments, c.census, c.summary,
			names.toArray(new String[names.size()]), references.toArray(new byte[references.size()][]),
			(header, model, bases)->{
				final String id=RrnaResourceIO.first(header), source=RrnaResourceIO.field(header, "source");
				final int genome=genomeIds.get(source);
				require(genome>=0 && !sourceIds.contains(id), "Every actual kept voter needs one original source and unique locus ID");
				final Genome owner=genomes.get(genome);
				final int next=groupSources.get(owner.group);
				sourceIds.put(id, sources.size());
				groupSources.put(owner.group, sources.size());
				sources.add(new Source(id, owner, bases.length, next, header));
				if(owner.group<0 || uncertainGroups.contains(owner.group)){unmappedDonor=true;}
				exact.add(bases);
			});
		require(sources.size()==checkedMembers.keptMembers(), "Only checked kept FASTA members may enter exact-sequence or source-locus evidence");
	}

	/** Corrected rows override historical raw bounds only for the same source ID. */
	private void positions(final String path, final int mode){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, mode==0 ? RAW_HEADER : mode==1 ? RECOVERED_HEADER : OVERLAY_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, mode==0 ? 21 : mode==1 ? 38 : 33);
				final int index=sourceIds.get(text(0));
				if(index<0){continue;}
				final Source source=sources.get(index);
				final String seqid=text(6);
				final int strand=strand(7);
				require(source.genome.id.equals(text(1)), "A corrected locus must retain its exact original donor genome");
				if(mode==0){
					require(!source.raw, "Duplicate original source-locus row");
					source.raw=true;
					final String intervals=text(8);
					final int dash=intervals.indexOf('-');
					if(dash>0 && intervals.indexOf(',', 0)<0 && intervals.indexOf(';', 0)<0 && intervals.indexOf('-', dash+1)<0){
						source.start=Integer.parseInt(intervals.substring(0, dash));
						source.stop=Integer.parseInt(intervals.substring(dash+1));
					}
				}else{
					require(!source.corrected && (mode!=1 || source.raw), "Corrected source coordinates must override exactly one appropriate original row");
					if(mode==2){require(p.termEquals("true", 23), "A kept clan member requires an admitted current overlay row");}
					final int offset=mode==1 ? 35 : 27;
					source.start=p.parseInt(offset);
					source.stop=p.parseInt(offset+1);
					require(p.parseInt(offset+2)==source.length && source.stop-source.start+1==source.length,
						"Final corrected source coordinates must describe the actual kept member length");
					require(source.headerSeqid==null || source.headerAgrees(seqid, strand),
						"Final source metadata must agree with retained FASTA header coordinates");
					source.corrected=true;
				}
				require(source.start>=0 && source.stop>=source.start && (source.seqid==null
					|| source.seqid.equals(seqid) && source.strand==strand), "Coordinate correction must preserve its contig and strand");
				source.seqid=seqid;
				source.strand=strand;
				if(mode==0 && source.headerSeqid!=null && !source.headerAgrees(seqid, strand)){source.needsCorrection=true;}
				source.known=source.start>0 && source.stop-source.start+1==source.length
					&& (!source.needsCorrection || source.corrected);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Saved Queries         ----------------*/
	/*--------------------------------------------------------------*/

	private void parents(){
		try(Input in=new Input(c.input+".parents.tsv.gz")){
			RrnaResourceIO.header(in.reader, NcrnaSavedRankSummary.PARENT_HEADER, c.input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 14);
				if(!p.termEquals("true", 9)){continue;}
				final int slot=p.parseInt(0), start=p.parseInt(5), stop=p.parseInt(6);
				require(slot>=0 && slot<samples.length && start>0 && stop>=start, "Target parents need a bound sample and valid one-based inclusive span");
				final Sample sample=samples[slot];
				require(sample.id==p.parseInt(1) && sample.path.equals(text(2)), "Parent input ID and exact path must match the accepted sample census");
				final Parent parent=new Parent(slot, sample, text(3), text(4), start, stop, strand(7), row);
				parent.next=sample.parents;
				sample.parents=parent;
				parents.add(parent);
			}
		}
		require(parents.size()==c.targets, "Retain the complete declared target inventory, including unrankable loci");
	}

	private void common(){
		int count=0;
		try(Input in=new Input(c.input+".loci.tsv.gz")){
			RrnaResourceIO.header(in.reader, NcrnaSavedRankSummary.LOCUS_HEADER, c.input);
			for(final Parent parent:parents){
				for(int method=0; method<c.methods; method++){
					final byte[] row=in.reader.nextLine();
					require(row!=null, "Saved locus-method blocks cannot be truncated before overlap classification");
					parse(row, 12);
					require(p.parseInt(0)==parent.slot && p.parseInt(1)==parent.sample.id && matches(2, parent.idBytes)
						&& p.parseInt(3)==method, "Overlap strata must preserve the original parent and complete method order");
					if(method==0){
						final int rank=p.parseInt(8);
						require(rank==-1 || rank>0 && rank<=c.models, "Common loci require an actual positive native rank");
						parent.common=rank>0;
						if(parent.common){count++;}
					}
				}
			}
			require(in.reader.nextLine()==null && count==c.common, "The seen/unseen comparison must conserve the original common native-rankable set");
		}
	}

	private void windows(){
		final ArrayList<Window> windows=new ArrayList<Window>();
		final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class);
		try(Input in=new Input(c.input+".windows.tsv.gz")){
			RrnaResourceIO.header(in.reader, WINDOW_HEADER, c.input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 15);
				final int slot=p.parseInt(0), start=p.parseInt(10), stop=p.parseInt(11), strand=p.parseInt(6);
				require(slot>=0 && slot<samples.length && p.parseInt(1)==samples[slot].id && start>=0
					&& stop>=start && stop<p.parseInt(12) && stop-start+1==p.parseInt(14) && (strand==0 || strand==1),
					"Saved windows use complete zero-based INCLUSIVE forward bounds and a declared orientation");
				final String key=slot+"_"+p.parseLong(2);
				require(!ids.contains(key), "A saved event must identify exactly one captured FASTA record");
				ids.put(key, windows.size());
				windows.add(new Window(samples[slot], text(4), start, stop, strand));
			}
		}
		sequences(c.input+".windows.fa.gz", read->{
			final int index=ids.get(read.id);
			require(index>=0 && read.bases!=null, "Every captured FASTA record needs its exact saved event metadata");
			final Window window=windows.get(index);
			require(!window.seen && read.length()==window.stop-window.start+1, "Duplicate or truncated captured window cannot supply gene sequence");
			window.seen=true;
			for(Parent parent=window.sample.parents; parent!=null; parent=parent.next){
				if(!parent.seqid.equals(window.seqid) || parent.start-1<window.start || parent.stop-1>window.stop){continue;}
				parent.captured=true;
				final byte[] candidate=extract(read.bases, window.start, window.stop, window.strand, parent.start, parent.stop, parent.strand);
				if(candidate==null){continue;}
				require(parent.sequence==null || Arrays.equals(parent.sequence, candidate),
					"Conflicting complete captured sequences cannot determine a label by read order: "+parent.id);
				if(parent.sequence==null){parent.sequence=candidate;}
			}
		});
		for(final Window window:windows){require(window.seen, "Captured window metadata without its FASTA record cannot silently disappear");}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Classification         ----------------*/
	/*--------------------------------------------------------------*/

	private int sourceState(final Parent parent){
		if(parent.sample.group<0){return UNKNOWN;}
		boolean unknown=unmappedDonor || uncertainGroups.contains(parent.sample.group);
		for(int i=groupSources.get(parent.sample.group); i>=0; i=sources.get(i).next){
			final Source source=sources.get(i);
			if(!source.known){unknown=true; continue;}
			if(!source.seqid.equals(parent.seqid)){
				if(!source.genome.path.equals(parent.sample.path)){unknown=true;}
				continue;
			}
			if(sameLocus(parent.start, parent.stop, parent.strand, source.start, source.stop, source.strand)){
				parent.matched=source.id;
				return SEEN;
			}
		}
		return unknown ? UNKNOWN : NOT_SEEN;
	}

	private void write(){
		final int[][] counts=new int[3][8];
		try(Output detail=new Output(c.output+".loci.tsv"); Output labels=new Output(c.output+".labels.tsv"); Output summary=new Output(c.output+".summary.tsv")){
			detail.writer.println(NcrnaSavedRankSummary.PARENT_HEADER+"\tcommon_native_rankable\tresolved_assembly_group\tcomplete_defined_query\texact_sequence\tsame_source_locus\tseen_stratum\tmatched_source_locus");
			labels.writer.println(NcrnaSavedRankSummary.STRATA_HEADER);
			for(final Parent parent:parents){
				final int sequence=parent.sequence==null ? UNKNOWN : exact.contains(parent.sequence) ? SEEN : NOT_SEEN;
				final int locus=sourceState(parent), state=classify(sequence, locus);
				final int[] row=counts[state];
				row[0]++;
				if(parent.common){
					row[1]++;
					if(sequence==SEEN){row[2]++;}
					if(locus==SEEN){row[3]++;}
					if(sequence==SEEN && locus==SEEN){row[4]++;}
					if(sequence==UNKNOWN){row[5]++;}
					if(locus==UNKNOWN){row[6]++;}
					if(parent.sequence!=null){row[7]++;}
				}
				detail.writer.print(b.clear().append(parent.original).tab().append(parent.common).tab().append(parent.sample.group)
					.tab().append(parent.sequence!=null).tab().append(STATES[sequence]).tab().append(STATES[locus]).tab().append(STATES[state])
					.tab().append(parent.matched==null ? "." : parent.matched).nl());
				labels.writer.print(b.clear().append(parent.slot).tab().append(parent.sample.id).tab()
					.append(parent.idBytes).tab().append(STATES[state]).nl());
			}
			summary.writer.println("seen_stratum\ttarget_loci\tcommon_loci\tcommon_exact_sequence_seen\tcommon_source_locus_seen\tcommon_both_seen\tcommon_sequence_unknown\tcommon_source_locus_unknown\tcommon_complete_defined_sequence");
			for(int state=0; state<counts.length; state++){
				b.clear().append(STATES[state]);
				for(final int count:counts[state]){b.tab().append(count);}
				summary.writer.print(b.nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns normalized sense-oriented DNA only under full inclusive containment. */
	public static byte[] extract(final byte[] window, final int windowStart0, final int windowStop0, final int windowStrand,
			final int start1, final int stop1, final int strand){
		require(window!=null && windowStart0>=0 && windowStop0-windowStart0+1==window.length
			&& start1>0 && stop1>=start1 && (windowStrand==0 || windowStrand==1) && (strand==0 || strand==1),
			"Extraction requires intact inclusive window coordinates and explicit query/window strands");
		if(start1-1<windowStart0 || stop1-1>windowStop0){return null;}
		final int length=stop1-start1+1, offset=windowStrand==0 ? start1-1-windowStart0 : windowStop0-(stop1-1);
		final byte[] result=new byte[length];
		for(int i=0; i<length; i++){
			final int base=RrnaResourceIO.base(window[offset+i]);
			if(base<0){return null;}
			result[windowStrand==strand ? i : length-1-i]=DNA[windowStrand==strand ? base : 3-base];
		}
		return result;
	}

	/** Presence relation, not one-to-one grading: same strand and reciprocal >half. */
	public static boolean sameLocus(final int start, final int stop, final int strand,
			final int otherStart, final int otherStop, final int otherStrand){
		require(start>0 && stop>=start && otherStart>0 && otherStop>=otherStart
			&& (strand==0 || strand==1) && (otherStrand==0 || otherStrand==1), "Source-locus overlap requires valid one-based inclusive spans and explicit strands");
		final long overlap=Math.min(stop, otherStop)-Math.max(start, otherStart)+1L;
		return strand==otherStrand && overlap>0 && 2*overlap>stop-start+1L && 2*overlap>otherStop-otherStart+1L;
	}

	public static int classify(final int exact, final int source){
		require(exact>=0 && exact<=UNKNOWN && source>=0 && source<=UNKNOWN, "Overlap components must be proven seen, proven absent or unknown");
		return exact==SEEN || source==SEEN ? SEEN : exact==NOT_SEEN && source==NOT_SEEN ? NOT_SEEN : UNKNOWN;
	}

	/** Drain native batches after a consumer failure; native close is not cancellation. */
	private static void sequences(final String path, final Consumer<Read> consumer){
		Throwable failure=null;
		try(SequenceInput in=new SequenceInput(StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1))){
			in.stream.start();
			for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
				for(final Read read:batch.list){
					if(failure!=null){continue;}
					try{
						require(read.mate==null, "Saved construction/window inputs must be unpaired FASTA records");
						consumer.accept(read);
					}catch(RuntimeException | Error invalid){failure=invalid;}
				}
				in.stream.returnList(batch);
			}
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			if(failure instanceof Error){throw (Error)failure;}
		}
	}

	private void parse(final byte[] row, final int fields){
		p.set(row);
		require(p.terms()==fields, "Saved overlap input field count changed; expected="+fields);
	}
	private String text(final int field){
		p.setBounds(field);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 would change source/query identity");
		return value;
	}
	/** Compare a retained UTF-8 identity without allocating per method row. */
	private boolean matches(final int field, final byte[] expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			if(p.line()[p.a()+i]!=expected[i]){return false;}
		}
		return true;
	}
	private int strand(final int field){
		require(p.termEquals("+", field) || p.termEquals("-", field), "Source/query strand must be explicit");
		return p.termEquals("+", field) ? 0 : 1;
	}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String input, output, members, ledger, assignments;
		public String census, summary, references, genomes, groups, assemblies;
		public String proxy, recovered;
		public String[] overlays;
		public int inputs=1000, genomeCount=39979, targets=2156, common=2007, models=508, methods=148;
	}
	private static final class Genome{
		Genome(final String id_, final String path_, final int group_){
			id=id_;
			path=path_;
			group=group_;
		}
		final String id, path;
		final int group;
	}
	private static final class Sample{
		Sample(final int id_, final String path_, final int group_){
			id=id_;
			path=path_;
			group=group_;
		}
		final int id, group;
		final String path;
		Parent parents;
	}
	private static final class Source{
		Source(final String id_, final Genome genome_, final int length_, final int next_, final String header){
			id=id_;
			genome=genome_;
			length=length_;
			next=next_;
			headerSeqid=RrnaResourceIO.optionalField(header, "contig");
			final String start=RrnaResourceIO.optionalField(header, "start");
			final String stop=RrnaResourceIO.optionalField(header, "stop");
			final String strand=RrnaResourceIO.optionalField(header, "strand");
			final boolean any=headerSeqid!=null || start!=null || stop!=null || strand!=null;
			require(!any || headerSeqid!=null && start!=null && stop!=null && ("+".equals(strand) || "-".equals(strand)),
				"Retained final-coordinate witnesses must contain contig, start, stop and strand together");
			headerStart=start==null ? 0 : Integer.parseInt(start);
			headerStop=stop==null ? 0 : Integer.parseInt(stop);
			headerStrand="-".equals(strand) ? 1 : 0;
			needsCorrection=RrnaResourceIO.optionalField(header, "recovery")!=null
				|| RrnaResourceIO.optionalField(header, "coordinate_model")!=null
				|| id.startsWith("RF01482_") || id.startsWith("RF01689_");
			require(!needsCorrection || any, "Corrected construction members must retain their final header coordinates");
			require(!any || headerStart>0 && headerStop>=headerStart && headerStop-headerStart+1==length,
				"A final construction header must describe its actual kept sequence length");
		}
		boolean headerAgrees(final String contig, final int strand){
			assert(headerSeqid!=null) : "Final-coordinate comparison requires the retained header witness checked by its caller";
			return headerSeqid.equals(contig) && headerStrand==strand && headerStart==start && headerStop==stop;
		}
		final String id;
		final Genome genome;
		final int length, next;
		final String headerSeqid;
		final int headerStart, headerStop, headerStrand;
		String seqid;
		int start, stop, strand;
		boolean raw, corrected, known, needsCorrection;
	}
	private static final class Parent{
		Parent(final int slot_, final Sample sample_, final String id_, final String seqid_, final int start_, final int stop_, final int strand_, final byte[] original_){
			slot=slot_;
			sample=sample_;
			id=id_;
			idBytes=id_.getBytes(StandardCharsets.UTF_8);
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			original=original_;
		}
		final int slot, start, stop, strand;
		final Sample sample;
		final String id, seqid;
		final byte[] original, idBytes;
		Parent next;
		String matched;
		byte[] sequence;
		boolean common, captured;
	}
	private static final class Window{
		Window(final Sample sample_, final String seqid_, final int start_, final int stop_, final int strand_){
			sample=sample_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
		}
		final Sample sample;
		final String seqid;
		final int start, stop, strand;
		boolean seen;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private final Config c;
	private final Sample[] samples;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final ObjectIntMap<String> genomeIds=new ObjectIntMap<String>(String.class), sourceIds=new ObjectIntMap<String>(String.class);
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Parent> parents=new ArrayList<Parent>();
	private final IntHashMap groupSources=new IntHashMap();
	private final IntHashMap uncertainGroups=new IntHashMap();
	private final NcrnaExactSequenceSet exact=new NcrnaExactSequenceSet();
	private NcrnaDdlMembers checkedMembers;
	private boolean unmappedDonor;
	public static final int NOT_SEEN=0, SEEN=1, UNKNOWN=2;
	private static final String[] STATES={"NOT_SEEN", "SEEN", "UNKNOWN"};
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
	public static final String ASSEMBLY_HEADER="input_slot\tinput_id\tfasta\tresolved_assembly_group\tkept_member_contributions\tassembly_status\tinterpretation";
	public static final String WINDOW_HEADER="input_slot\tinput_id\tevent\tcontig_ordinal\tseqid\theader\tstrand0\tpass\toriented_start0\toriented_stop0\tforward_start0\tforward_stop0\tcontig_length\tsource\twindow_length";
	public static final String RAW_HEADER="locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines";
	public static final String RECOVERED_HEADER=RAW_HEADER+"\tcm_state\traw_hits\toverlapping_sense_hits\tcm_model_from\tcm_model_to\tcm_seq_from\tcm_seq_to\tcm_trunc\tcm_bits\toverlap_bp\tleft_excess\tright_excess\tgeometry\trecovery\tcore_start1\tcore_stop1\tcore_length";
	public static final String OVERLAY_HEADER=RAW_HEADER+"\tsource_family\ttruth_subtype\tadmitted\tsupport_models\tcoordinate_model\tevidence_stage\tcurrent_start1\tcurrent_stop1\tcore_length\tdiagnostic_length_band\tmember_family\tgeometry_status";
}
