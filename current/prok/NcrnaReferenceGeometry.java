package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import map.LongHashSet;
import map.LongIntMap;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;
import shared.Tools;

/** Descriptive connected groups from saved model hits and their saved overlap edges.
 * Groups are review containers, not genes: transitive membership does not imply
 * pairwise overlap. Invalid intervals remain ungrouped; opposite-strand edges
 * remain diagnostic. No score, span floor, subtype winner or admission is selected.
 * Original input/table/line keys preserve access to native evidence, including
 * E-values, which are not duplicated here. One thread owns all retained state.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceGeometry{

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	static void process(final String input, final String output, final NcrnaEvidenceCatalog catalog, final Roster roster){
		load(input, catalog, roster, true).write(output);
	}

	/** Shared native event/geometry loader for downstream descriptive evidence joins. */
	static NcrnaReferenceGeometry load(final String input, final NcrnaEvidenceCatalog catalog, final Roster roster){
		return load(input, catalog, roster, false);
	}
	/** Explicit opt-in protects legacy decision tools from silently changing GA semantics. */
	static NcrnaReferenceGeometry load(final String input, final NcrnaEvidenceCatalog catalog, final Roster roster, final boolean allowBatch){
		final NcrnaReferenceProfiles profiles=NcrnaReferenceProfiles.load(input, catalog, roster, allowBatch);
		final NcrnaReferenceGeometry geometry=new NcrnaReferenceGeometry(catalog, roster, profiles);
		geometry.readHits(input+".hits.tsv.gz");
		geometry.readPairs(input+".overlaps.tsv.gz");
		geometry.group();
		return geometry;
	}

	Hit[] hits(){return hits;}
	Hit[] ordered(){return ordered;}
	ArrayList<Pair> pairs(){return pairs;}
	NcrnaReferenceProfiles profiles(){return profiles;}
	int groups(){return groups;}
	InputState input(final String id){return inputs.get(id);}

	private NcrnaReferenceGeometry(final NcrnaEvidenceCatalog catalog_, final Roster roster_, final NcrnaReferenceProfiles profiles_){
		catalog=catalog_;
		roster=roster_;
		profiles=profiles_;
		require(profiles.events()<=Integer.MAX_VALUE, "Retained event arrays require an int-sized event population");
		hits=new Hit[(int)profiles.events()];
		parent=new int[hits.length];
		rank=new byte[hits.length];
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			inputs.put(genome.id, new InputState(genome, i, profiles.selected(genome.id)));
		}
		families=new int[roster.queried()];
		final ObjectIntMap<String> names=new ObjectIntMap<String>(String.class);
		for(int i=0; i<families.length; i++){
			final String name=roster.query(i).scope.family;
			if(!names.contains(name)){names.put(name, names.size());}
			families[i]=names.get(name);
		}
		familyCount=names.size();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved event join       ----------------*/
	/*--------------------------------------------------------------*/

	private void readHits(final String path){
		final LineParser1 p=new LineParser1('\t');
		final ObjectMap<Target,Target> targets=new ObjectMap<Target,Target>(Target.class, Target.class);
		final Target probe=new Target();
		InputState input=null;
		Target target=null;
		int query=-1, count=0;
		final long[] modelCounts=new long[roster.queried()];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, profiles.batchEvidence() ? BuildNcrnaReference.BATCH_HIT_HEADER : BuildNcrnaReference.HIT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==(profiles.batchEvidence() ? 24 : 21) && count<hits.length, "Saved hit schema/count differs from its summary: "+path);
				if(input==null || !p.termEquals(input.genome.id, 0)){input=inputs.get(p.parseString(0));}
				require(input!=null && input.selected, "Saved hit names an input outside the observed catalog population");
				final Report report=input.report(p, 1, true);
				require(profiles.batchEvidence() || input.reports.size()==1, "One input must retain its original native table identity");
				if(target==null || !same(p, 3, target.bytes)){
					probe.set(p, 3);
					target=targets.get(probe);
					if(target==null){
						target=probe.copy();
						targets.put(target, target);
					}
				}
				if(query<0 || !p.termEquals(roster.query(query).accession, 7)){query=findModel(p);}
				final Model model=roster.query(query);
				require(p.termEquals(model.release, 8) && p.termEquals(model.scope.family, 9)
					&& p.termEquals(model.scope.role, 10) && p.parseInt(13)==model.length, "Saved hit model metadata differs from the supplied roster");
				final Hit hit=new Hit(input, target, report.path.bytes, query, count, p, model, profiles.batchEvidence());
				require(!report.lines.contains(hit.line), "Repeated input/table/native-line event");
				report.lines.put(hit.line, count);
				hits[count]=hit;
				parent[count]=count;
				count++;
				modelCounts[query]++;
			}
		}
		require(count==hits.length, "Saved events must conserve the original summary total");
		for(int i=0; i<modelCounts.length; i++){
			require(modelCounts[i]==profiles.modelEvents(i), "Saved event count differs from its accession histogram: "+roster.query(i).accession);
		}
	}

	/** Small fixed roster search on model changes; no transient field Strings. */
	private int findModel(final LineParser1 p){
		for(int i=0; i<roster.queried(); i++){if(p.termEquals(roster.query(i).accession, 7)){return i;}}
		throw new IllegalArgumentException("Saved hit uses an unknown or unqueried model");
	}

	private void readPairs(final String path){
		final LineParser1 p=new LineParser1('\t');
		final LongHashSet seen=new LongHashSet(16, 0.7f);
		InputState input=null;
		long historical=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.OVERLAP_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==9, "Saved overlap rows require nine fields");
				if(p.termEquals("HISTORICAL_MODEL", 1)){
					historical++;
					continue;
				}
				require(p.termEquals("MODEL_MODEL", 1), "Unknown saved overlap kind");
				if(input==null || !p.termEquals(input.genome.id, 0)){input=inputs.get(p.parseString(0));}
				require(input!=null && input.selected, "Model pair needs an observed input with retained hits");
				final int ai=input.find(p, 2, p.parseLong(3)), bi=input.find(p, 4, p.parseLong(5));
				require(ai>=0 && bi>=0 && ai!=bi, "Model pair must name two distinct retained events");
				final long key=((long)Math.min(ai, bi)<<32)|(Math.max(ai, bi)&0xffffffffL);
				require(seen.add(key), "Repeated saved model pair");
				final Hit a=hits[ai], b=hits[bi];
				final long overlap=overlap(a, b);
				require(a.target==b.target && overlap>0 && overlap==p.parseLong(6), "Saved pair target/overlap disagrees with its event coordinates");
				final boolean sameStrand=a.strand==b.strand;
				require(p.termEquals(sameStrand ? "SAME_STRAND_GEOMETRY" : "OPPOSITE_STRAND_GEOMETRY", 7)
					&& p.termEquals("false", 8), "Model geometry cannot carry a historical support-winner flag");
				pairs.add(new Pair(a, b));
				if(!a.valid() || !b.valid()){
					invalidPairs++;
				}else if(sameStrand){
					samePairs++;
					join(ai, bi);
					a.neighbor(bi);
					b.neighbor(ai);
				}else{oppositePairs++;}
			}
		}
		require(pairs.size()==profiles.modelPairs() && historical==profiles.historicalPairs(), "Saved overlap counts disagree with the completed evidence summary");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Connected groups       ----------------*/
	/*--------------------------------------------------------------*/

	/** Union-find uses only actual, valid same-strand edges; it selects no representative. */
	private int root(int index){
		assert(index>=0 && index<parent.length) : "Pair lookup must resolve a retained event before union-find";
		while(parent[index]!=index){
			parent[index]=parent[parent[index]];
			index=parent[index];
		}
		return index;
	}
	private void join(final int left, final int right){
		final int a=root(left), b=root(right);
		if(a==b){return;}
		if(rank[a]<rank[b]){parent[a]=b;}
		else{
			parent[b]=a;
			if(rank[a]==rank[b]){rank[a]++;}
		}
	}
	private void group(){
		ordered=hits.clone();
		Arrays.sort(ordered, POSITION);
		final int[] numbers=new int[hits.length];
		for(final Hit hit:ordered){
			if(!hit.valid()){
				invalidHits++;
				continue;
			}
			final int root=root(hit.index);
			if(numbers[root]==0){numbers[root]=++groups;}
			hit.group=numbers[root];
		}
		// Group order also makes output robust to input row order and union-find roots.
		Arrays.sort(ordered, GROUP_ORDER);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reports            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		writeGroups(prefix+".groups.tsv.gz");
		writeMembers(prefix+".group_members.tsv.gz");
		writePairs(prefix+".model_pairs.tsv.gz");
		writeBridges(prefix+".bridges.tsv.gz");
		try(Output out=new Output(prefix+".geometry_summary.tsv")){
			out.writer.println("inputs_observed\tcatalog_inputs\tqueried_models\traw_hit_events\tgrouped_events\tunassessable_events\treview_groups\tmodel_pairs\tvalid_same_strand_pairs\tvalid_opposite_strand_pairs\tunassessable_pairs\tbridge_witness_rows\tevidence_status\tstatus");
			final ByteBuilder b=new ByteBuilder();
			b.append(profiles.observed()).tab().append(catalog.size()).tab().append(roster.queried()).tab().append(hits.length).tab();
			b.append(hits.length-invalidHits).tab().append(invalidHits).tab().append(groups).tab().append(pairs.size()).tab();
			b.append(samePairs).tab().append(oppositePairs).tab().append(invalidPairs).tab().append(bridges).tab().append(profiles.evidenceState()).tab();
			out.writer.print(b.append("DESCRIPTIVE_GEOMETRY_NOT_GENES").nl());
		}
	}

	private void writeGroups(final String path){
		final int[] seenModel=new int[roster.queried()], seenFamily=new int[familyCount];
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println("group_id\tinput_id\ttarget\tstrand\tstart\tstop\tevents\tfull_model_untruncated\tshorter_untruncated\tnative_truncated\tleft_edge_events\tright_edge_events\tdistinct_models\tdistinct_families\tstatus");
			int i=0;
			while(i<ordered.length){
				final Hit first=ordered[i];
				if(first.group==0){
					i++;
					continue;
				}
				int end=i, stop=first.stop, full=0, shorter=0, truncated=0, left=0, right=0, models=0, families=0;
				while(end<ordered.length && ordered[end].group==first.group){
					final Hit hit=ordered[end++];
					assert(hit.input==first.input && hit.target==first.target && hit.strand==first.strand) : "Only valid same-context pair edges may form a review group";
					stop=Math.max(stop, hit.stop);
					if(hit.full()){full++;}
					else if(hit.truncation==0){shorter++;}
					else{truncated++;}
					if(hit.left){left++;}
					if(hit.right){right++;}
					if(seenModel[hit.query]!=hit.group){
						seenModel[hit.query]=hit.group;
						models++;
					}
					final int family=this.families[hit.query];
					if(seenFamily[family]!=hit.group){
						seenFamily[family]=hit.group;
						families++;
					}
				}
				b.clear().append(first.group).tab().append(first.input.genome.id).tab().append(first.target.bytes).tab().append(first.strand).tab();
				b.append(first.start).tab().append(stop).tab().append(end-i).tab().append(full).tab().append(shorter).tab().append(truncated).tab();
				out.writer.print(b.append(left).tab().append(right).tab().append(models).tab().append(families).append("\tOVERLAP_COMPONENT_NOT_GENE\n"));
				i=end;
			}
		}
	}

	private void writeMembers(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println("group_id\tinput_id\tnative_table\tnative_line\ttarget\tstart\tstop\tstrand\taccession\treference_family\tevidence_role\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tbits\tGA_pass\trange\tleft_contig_edge\tright_contig_edge"
				+(profiles.batchEvidence() ? "\tnative_included\tGA_decision_basis" : ""));
			for(final Hit hit:ordered){
				final Model model=roster.query(hit.query);
				groupId(b.clear(), hit).tab().append(hit.input.genome.id).tab().append(hit.table).tab().append(hit.line).tab();
				b.append(hit.target.bytes).tab().append(hit.start).tab().append(hit.stop).tab().append(hit.strand).tab().append(model.accession).tab();
				b.append(model.scope.family).tab().append(model.scope.role).tab().append(hit.from).tab().append(hit.to).tab().append(model.length).tab();
				b.append(hit.truncation).tab().append(hit.bits, 6).tab().append(hit.ga).tab().append(RANGES[hit.range]).tab();
				b.append(hit.left).tab().append(hit.right);
				if(profiles.batchEvidence()){b.tab().append(hit.included).tab().append(BuildNcrnaReference.BOUND_GA_BASIS);}
				out.writer.print(b.nl());
			}
		}
	}

	private void writePairs(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println((profiles.batchEvidence() ? "input_id\tleft_native_table\tleft_line\tright_native_table\tright_line" : "input_id\tnative_table\tleft_line\tright_line")
				+"\tleft_accession\tright_accession\tleft_group\tright_group\trelation\toverlap_bases\tleft_length\tright_length\tleft_overlap_fraction\tright_overlap_fraction\tcontext\tboth_full_model_untruncated");
			for(final Pair pair:pairs){
				final Hit a=pair.a, z=pair.b;
				final long overlap=overlap(a, z);
				b.clear().append(a.input.genome.id).tab().append(a.table).tab().append(a.line).tab();
				if(profiles.batchEvidence()){b.append(z.table).tab();}
				b.append(z.line).tab();
				b.append(roster.query(a.query).accession).tab().append(roster.query(z.query).accession).tab();
				groupId(b, a).tab();
				groupId(b, z).tab().append(relation(a, z)).tab().append(overlap).tab().append(a.length()).tab().append(z.length()).tab();
				b.append(overlap/(double)a.length(), 6).tab().append(overlap/(double)z.length(), 6).tab();
				b.append(!a.valid() || !z.valid() ? "UNASSESSABLE" : a.strand==z.strand ? "SAME_STRAND" : "OPPOSITE_STRAND").tab();
				out.writer.print(b.append(a.full() && z.full()).nl());
			}
		}
	}

	/** One deterministic disjoint pair per bridge and witness family, never a gene call.
	 * All pair edges and members survive separately; this is a compact witness, not
	 * an exhaustive triple enumeration. Scratch arrays are reused across bridges.
	 */
	private void writeBridges(final String path){
		Hit[] buffer=new Hit[16];
		final Hit[] earliest=new Hit[familyCount];
		final boolean[] emitted=new boolean[familyCount];
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println((profiles.batchEvidence() ? "group_id\tinput_id\tbridge_native_table\tbridge_line\tleft_native_table\tleft_line\tright_native_table\tright_line"
				: "group_id\tinput_id\tnative_table\tbridge_line\tleft_line\tright_line")
				+"\tbridge_accession\tleft_accession\tright_accession\twitness_family\tbridge_left_overlap\tbridge_right_overlap\tdisjoint_gap_bases\tcovers_outer_span\tstatus");
			for(final Hit bridge:ordered){
				if(bridge.neighbors==null || bridge.neighbors.size<2){continue;}
				if(buffer.length<bridge.neighbors.size){buffer=new Hit[bridge.neighbors.size];}
				int n=0;
				for(int i=0; i<bridge.neighbors.size; i++){
					final Hit hit=hits[bridge.neighbors.get(i)];
					if(hit.full()){buffer[n++]=hit;}
				}
				Arrays.sort(buffer, 0, n, POSITION);
				Arrays.fill(earliest, null);
				Arrays.fill(emitted, false);
				for(int i=0; i<n; i++){
					final Hit z=buffer[i];
					final int family=families[z.query];
					final Hit a=earliest[family];
					if(a!=null && a.stop<z.start && !emitted[family]){
						assert(overlap(bridge, a)>0 && overlap(bridge, z)>0) : "Witnesses must be retained same-strand neighbors of the bridge";
						b.clear().append(bridge.group).tab().append(bridge.input.genome.id).tab().append(bridge.table).tab().append(bridge.line).tab();
						if(profiles.batchEvidence()){b.append(a.table).tab();}
						b.append(a.line).tab();
						if(profiles.batchEvidence()){b.append(z.table).tab();}
						b.append(z.line).tab().append(roster.query(bridge.query).accession).tab();
						b.append(roster.query(a.query).accession).tab().append(roster.query(z.query).accession).tab().append(roster.query(z.query).scope.family).tab();
						b.append(overlap(bridge, a)).tab().append(overlap(bridge, z)).tab().append((long)z.start-a.stop-1).tab();
						out.writer.print(b.append(bridge.start<=a.start && bridge.stop>=z.stop).append("\tCANDIDATE_BRIDGE_NOT_FUSION\n"));
						emitted[family]=true;
						bridges++;
					}
					if(a==null || z.stop<a.stop){earliest[family]=z;}
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Retained state       ----------------*/
	/*--------------------------------------------------------------*/

	static final class InputState{
		InputState(final Genome genome_, final int index_, final boolean selected_){
			genome=genome_;
			index=index_;
			selected=selected_;
		}
		final Genome genome;
		final int index;
		final boolean selected;
		/** Resolves a literal encoded report path plus physical line without field Strings. */
		int find(final LineParser1 p, final int tableField, final long line){
			final Report report=report(p, tableField, false);
			return report==null ? -1 : report.lines.get(line);
		}
		private Report report(final LineParser1 p, final int field, final boolean create){
			lookup.set(p, field);
			Report report=reports.get(lookup);
			if(report==null && create){
				require(lookup.to>lookup.from, "Native report identity cannot be empty");
				report=new Report(lookup.copy());
				reports.put(report.path, report);
			}
			return report;
		}
		private final ObjectMap<Target,Report> reports=new ObjectMap<Target,Report>(Target.class, Report.class);
		private final Target lookup=new Target();
	}
	/** One retained report identity owns its native-line index within an input. */
	private static final class Report{
		Report(final Target path_){path=path_;}
		final Target path;
		final LongIntMap lines=new LongIntMap(16);
	}

	/** Byte-key interning uses one mutable lookup probe; only retained keys are copied.
	 * Stored targets are never mutated. Equal bytes in different inputs may share
	 * this name object but do not share genomic identity or connected groups.
	 */
	static final class Target{
		void set(final LineParser1 p, final int field){
			p.setBounds(field);
			set(p.line(), p.a(), p.b());
		}
		/** Borrows a byte range for lookup; only copy() may become a stored key. */
		void set(final byte[] data, final int start, final int end){
			assert(start>=0 && end>=start && end<=data.length) : "Target-key hashing must stay within the supplied byte range";
			bytes=data;
			from=start;
			to=end;
			hash=1;
			for(int i=from; i<to; i++){hash=31*hash+bytes[i];}
		}
		Target copy(){
			final Target result=new Target();
			result.bytes=Arrays.copyOfRange(bytes, from, to);
			result.to=result.bytes.length;
			result.hash=hash;
			return result;
		}
		@Override
		public int hashCode(){return hash;}
		@Override
		public boolean equals(final Object object){
			if(!(object instanceof Target)){return false;}
			final Target other=(Target)object;
			if(to-from!=other.to-other.from){return false;}
			for(int i=0; i<to-from; i++){if(bytes[from+i]!=other.bytes[other.from+i]){return false;}}
			return true;
		}
		byte[] bytes;
		int from, to, hash;
	}

	static final class Hit{
		Hit(final InputState input_, final Target target_, final byte[] table_, final int query_, final int index_, final LineParser1 p, final Model model, final boolean bound){
			input=input_;
			target=target_;
			table=table_;
			query=query_;
			index=index_;
			line=p.parseLong(2);
			start=p.parseInt(4);
			stop=p.parseInt(5);
			strand=p.termEquals("+", 6) ? '+' : p.termEquals("-", 6) ? '-' : '?';
			from=p.parseInt(11);
			to=p.parseInt(12);
			truncation=p.parseInt(14);
			bits=p.parseDouble(15);
			ga=bool(p, 17);
			included=bound && bool(p, 22);
			require(!bound || ga==included && p.termEquals(BuildNcrnaReference.BOUND_GA_BASIS, 23), "Batch GA membership must retain its verified native inclusion basis");
			int state=-1;
			for(int i=0; i<RANGES.length; i++){if(p.termEquals(RANGES[i], 18)){state=i;}}
			range=state;
			left=bool(p, 19);
			right=bool(p, 20);
			require(line>0 && start>0 && stop>=start && strand!='?' && from>0 && to>=from && to<=model.length
				&& truncation>=0 && truncation<=3 && Double.isFinite(bits) && range>=0, "Invalid saved hit coordinates/strand/model/score/range");
			require(range==0 || (!left && !right), "Unassessable intervals cannot have physical contig-edge evidence");
			full=from==1 && to==model.length && truncation==0;
		}
		boolean valid(){return range==0;}
		boolean full(){return full;}
		long length(){return (long)stop-start+1;}
		void neighbor(final int other){
			if(neighbors==null){neighbors=new IntList(2);}
			neighbors.add(other);
		}
		final InputState input;
		final Target target;
		final byte[] table;
		final int query, index, start, stop, from, to, truncation, range;
		final long line;
		final char strand;
		final double bits;
		/** included is meaningful only in an explicitly bound batch generation. */
		final boolean ga, included, left, right, full;
		int group;
		IntList neighbors;
	}
	static final class Pair{
		Pair(final Hit a_, final Hit b_){
			a=a_;
			b=b_;
		}
		final Hit a, b;
	}
	private static boolean same(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		if(expected==null || expected.length!=p.b()-p.a()){return false;}
		for(int i=0; i<expected.length; i++){if(expected[i]!=p.line()[p.a()+i]){return false;}}
		return true;
	}
	private static boolean bool(final LineParser1 p, final int field){
		if(p.termEquals("true", field)){return true;}
		require(p.termEquals("false", field), "Saved boolean must be true or false");
		return false;
	}
	private static long overlap(final Hit a, final Hit b){return Math.max(0L, (long)Math.min(a.stop, b.stop)-Math.max(a.start, b.start)+1);}
	private static String relation(final Hit a, final Hit b){
		if(a.start==b.start && a.stop==b.stop){return "EQUAL";}
		if(a.start<=b.start && a.stop>=b.stop){return "LEFT_CONTAINS_RIGHT";}
		if(b.start<=a.start && b.stop>=a.stop){return "RIGHT_CONTAINS_LEFT";}
		return "PARTIAL_OVERLAP";
	}
	private static ByteBuilder groupId(final ByteBuilder b, final Hit hit){return hit.group==0 ? b.append("NA") : b.append(hit.group);}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static final Comparator<Hit> POSITION=(a, b)->{
		int c=Integer.compare(a.input.index, b.input.index);
		if(c==0){c=Tools.compare(a.target.bytes, b.target.bytes);}
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		if(c==0){c=Tools.compare(a.table, b.table);}
		return c==0 ? Long.compare(a.line, b.line) : c;
	};
	private static final Comparator<Hit> GROUP_ORDER=(a, b)->a.group==b.group ? POSITION.compare(a, b) : Integer.compare(a.group, b.group);
	private final NcrnaEvidenceCatalog catalog;
	private final Roster roster;
	private final NcrnaReferenceProfiles profiles;
	private final ObjectMap<String,InputState> inputs=new ObjectMap<String,InputState>(String.class, InputState.class);
	private final ArrayList<Pair> pairs=new ArrayList<Pair>();
	private final Hit[] hits;
	private Hit[] ordered;
	private final int[] parent, families;
	private final byte[] rank;
	private final int familyCount;
	private int groups, invalidHits;
	private long samePairs, oppositePairs, invalidPairs, bridges;
	private static final String[] RANGES={"VALID_INTERVAL", "OUT_OF_RANGE", "UNRESOLVED_LITERAL_TARGET"};
	static final String[] OUTPUT_SUFFIXES={".groups.tsv.gz", ".group_members.tsv.gz", ".model_pairs.tsv.gz", ".bridges.tsv.gz", ".geometry_summary.tsv"};
}
