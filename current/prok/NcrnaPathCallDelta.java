package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectSet;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Source-qualified comparison of saved ordinary caller paths, one input at a time.
 * Exact feature multiplicity, coordinate links and physical reference overlap are
 * separate observations. Ambiguous overlaps never become guessed gene identities.
 * Computes physical reference-overlap coverage; no caller, sequence reader,
 * one-to-one reference grading or biological accuracy adjudication runs here.
 * @author Ganyu
 */
public final class NcrnaPathCallDelta{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPathCallDelta.class, false);
		try{
			new NcrnaPathCallDelta(pp.args).process();
			pp.outstream.println("NCRNA_PATH_CALL_DELTA_COMPLETE_SAVED_GFF_OBSERVATIONS_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaPathCallDelta(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("beforeconfig")){beforeConfig=b;}
			else if(a.equals("afterconfig")){afterConfig=b;}
			else if(a.equals("sources")){sources=b;}
			else if(a.equals("genomes")){genomes=b;}
			else if(a.equals("reference")){reference=b;}
			else if(a.equals("legacyfilenames")){legacyNames=Parse.parseBoolean(b);}
			else if(a.equals("expectedinputs")){expected=NcrnaNumbers.intKMG(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown path comparison option: "+arg);}
		}
		before=parser.in1;
		after=parser.in2;
		out=parser.out1;
		require(before!=null && after!=null && beforeConfig!=null && afterConfig!=null
			&& genomes!=null && reference!=null && out!=null && expected>0 && (sources!=null)!=legacyNames,
			"Bind both saved GFFs/configs, genomes=normalized_OFF_genomes, reference=normalized_OFF_references, expectedinputs=, sources= OR legacyfilenames=t, out=fresh_directory");
		require(!parser.overwrite && !parser.append && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Path comparison retains complete populations and never overwrites evidence");
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){
		bindings();
		references();
		require(new File(out).mkdir(), "Path comparison needs a fresh output directory: "+out);
		try(Cursor a=new Cursor(before);
			Cursor b=new Cursor(after);
			Output changes=output("changes.tsv.gz");
			Output links=output("overlap_links.tsv.gz");
			Output neighbours=output("neighbours.tsv.gz");
			Output coverage=output("reference_coverage.tsv");
			Output refLinks=output("changed_reference_links.tsv.gz");
			Output overlaps=output("clan_feature_overlaps.tsv.gz");
			Output bindings=output("bindings.tsv")){
			final Reports reports=new Reports(changes, links, neighbours, coverage, refLinks, overlaps);
			reports.headers();
			bindings.writer.println("input_ordinal0\tinput_id\tfasta");
			for(int ordinal=0; ordinal<order.size(); ordinal++){
				final Binding binding=order.get(ordinal);
				bindings.writer.print(new ByteBuilder().append(ordinal).tab().append(binding.id).tab().append(binding.fasta).nl());
				final ArrayList<Feature> left=a.input(ordinal);
				final ArrayList<Feature> right=b.input(ordinal);
				compare(binding, left, right, reports);
			}
			require(a.pending==null && b.pending==null, "Saved GFF contains an undeclared extra input");
		}
		try(Output summary=output("summary.tsv");
			Output refs=output("reference_summary.tsv");
			Output overlaps=output("clan_overlap_summary.tsv")){
			summary.writer.println("cohort\tbefore_features\tafter_features\texact_features\tremoved_exact_rows\tadded_exact_rows\tcoordinate_pairs_from\tcoordinate_pairs_to\tmetadata_pairs_from\tmetadata_pairs_to\tremoved_without_any_counterpart_overlap\tadded_without_any_counterpart_overlap\toverlapping_unpaired_removed\toverlapping_unpaired_added\texact_duplicate_extra_before\texact_duplicate_extra_after\tspatial_duplicate_extra_before\tspatial_duplicate_extra_after\tsplit_candidates_before\tmerge_candidates_after");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<COHORTS.length; i++){
				final long[] c=counts[i];
				require(c[0]==c[2]+c[3] && c[1]==c[2]+c[4], "Exact multiset accounting must conserve every feature: "+COHORTS[i]);
				require(c[3]==c[5]+c[7]+c[9]+c[11] && c[4]==c[6]+c[8]+c[10]+c[12],
					"Every changed row must have exactly one explanation: unique coordinate pair, unique metadata pair, no counterpart overlap, or unpaired overlap: "+COHORTS[i]);
				row.clear().append(COHORTS[i]);
				for(final long n:c){row.tab().append(n);}
				summary.writer.print(row.nl());
			}
			refs.writer.println("reference_scope\ttotal\tcovered_before\tcovered_after\tnewly_covered\tnewly_uncovered\tqualification");
			for(int i=0; i<referenceCounts.length; i++){
				final long[] c=referenceCounts[i];
				require(c[1]+c[3]-c[4]==c[2], "Reference coverage changes must conserve the same frozen reference IDs");
				row.clear().append(i==0 ? "ALL" : i==1 ? "SCORED" : "MASKED");
				for(final long n:c){row.tab().append(n);}
				refs.writer.print(row.append("\tSAME_SOURCE_CONTIG_STRAND_POSITIVE_OVERLAP_NOT_GRADED_ACCURACY\n"));
			}
			overlaps.writer.println("side\tclan_calls\taffected_clan_calls\tpositive_feature_pairs\tsummed_overlap_bp\tmaximum_overlap_bp\tmean_overlap_bp\tqualification");
			for(int side=0; side<2; side++){
				final long[] c=overlapCounts[side];
				require(c[1]<=counts[0][side] && (c[0]>0)==(c[1]>0), "Affected clan calls are a subset of the whole path and require a positive overlap");
				row.clear().append(side==0 ? "BEFORE" : "AFTER").tab().append(counts[0][side])
					.tab().append(c[1]).tab().append(c[0]).tab().append(c[2]).tab();
				if(c[0]>0){row.append(c[3]);}
				row.tab();
				if(c[0]>0){row.append(c[2]/(double)c[0], 6);}
				overlaps.writer.print(row.append("\tUNIQUE_UNORDERED_PHYSICAL_FEATURE_PAIRS_AT_LEAST_ONE_CLAN_CALL_EITHER_STRAND\n"));
			}
		}
		try(Output done=output("complete.tsv")){
			done.writer.println("inputs\tqualification");
			done.writer.print(new ByteBuilder().append(order.size()).append("\tALL_INPUTS_AND_FEATURE_MULTIPLICITIES_ACCOUNTED_NO_COUNT_CHANGE_STOP\n"));
		}
	}

	/** Configs independently bind every ordinal to its exact FASTA and saved output. */
	private void bindings(){
		final ArrayList<String> left=config(beforeConfig, before);
		final ArrayList<String> right=config(afterConfig, after);
		require(left.equals(right) && left.size()==expected, "Before/after configurations must have identical complete ordered input populations");
		final ObjectSet<String> wanted=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(genomes)){
			RrnaResourceIO.header(input.reader, "input_id\tcalls\tmatched\tmasked\tunmatched", genomes);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5 && p.length(0)>0 && wanted.add(p.parseString(0)), "Normalized OFF genomes must identify every input once");
			}
		}
		require(wanted.size()==expected, "Normalized OFF population differs from the declared input count");
		if(legacyNames){
			for(final String fasta:left){
				final String name=new File(fasta).getName();
				require(name.endsWith(".fna.gz"), "Explicit historical filename binding requires the recorded .fna.gz convention");
				addBinding(name.substring(0, name.length()-7), fasta, wanted);
			}
		}else{
			try(Input input=new Input(sources)){
				RrnaResourceIO.header(input.reader, "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff", sources);
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==5 && order.size()<expected && integer(p, 0)==order.size()
						&& p.termEquals(left.get(order.size()), 3), "Source roster must match the exact caller configuration order");
					addBinding(p.parseString(2), p.parseString(3), wanted);
				}
			}
		}
		require(order.size()==expected && byId.size()==wanted.size(), "Source IDs and OFF population must be a complete bijection");
	}
	private void addBinding(final String id, final String fasta, final ObjectSet<String> wanted){
		require(wanted.contains(id) && !byId.contains(id), "Unbound or duplicate source input ID: "+id);
		final Binding binding=new Binding(order.size(), id, fasta);
		byId.put(id, binding);
		order.add(binding);
	}
	private static ArrayList<String> config(final String path, final String gff){
		final ArrayList<String> result=new ArrayList<String>();
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		int outputs=0;
		boolean ordered=false;
		boolean qualified=false;
		try(ConfigInput input=new ConfigInput(path)){
			while(input.next()){
				final String key=input.key.toLowerCase(Locale.ROOT);
				final String value=input.value;
				if(key.equals("in") || key.equals("in1")){
					require(new File(value).isAbsolute() && seen.add(value), "Each input must have a unique absolute FASTA path");
					result.add(value);
				}else if(key.equals("out") || key.equals("out1")){
					require(++outputs==1 && value.equals(gff), "Config output must name this exact saved GFF");
				}else if(key.equals("ordered")){ordered=Parse.parseBoolean(value);}
				else if(key.equals("gffinputordinal")){qualified=Parse.parseBoolean(value);}
				else{require(!key.equals("config") && !key.equals("in2"), "Flat unpaired caller configurations are required");}
			}
		}
		require(outputs==1 && ordered && qualified, "Source-qualified comparison requires ordered=t and gffinputordinal=t");
		return result;
	}
	private void references(){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(reference)){
			RrnaResourceIO.header(input.reader, "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored", reference);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==7, "Reference evidence requires all seven normalized fields");
				final Binding binding=byId.get(p.parseString(0));
				require(binding!=null && seen.add(p.parseString(0)+'\t'+p.parseString(1)), "Every frozen reference belongs to one declared source and unique ID");
				binding.references.add(new Reference(p));
			}
		}
	}

	private void compare(final Binding binding, final ArrayList<Feature> a, final ArrayList<Feature> b, final Reports reports){
		Collections.sort(a);
		Collections.sort(b);
		final ObjectMap<String, Group> ac=contigs(a);
		final ObjectMap<String, Group> bc=contigs(b);
		inventory(a, ac, 0);
		inventory(b, bc, 1);
		int i=0;
		int j=0;
		while(i<a.size() && j<b.size()){
			final Feature x=a.get(i);
			final Feature y=b.get(j);
			final int c=x.compareTo(y);
			if(c==0){
				x.exact=true;
				y.exact=true;
				add(x, 2);
				i++;
				j++;
			}else if(c<0){i++;}
			else{j++;}
		}
		pathOverlaps(binding, a, ac, 0, reports);
		pathOverlaps(binding, b, bc, 1, reports);
		for(final Feature x:a){
			if(x.exact){continue;}
			final Group group=group(bc, x.seqid);
			for(int k=group.first(x.start); k<group.rows.size(); k++){
				final Feature y=group.rows.get(k);
				if(y.start>x.stop){break;}
				if(!y.exact && compatible(x, y) && reciprocal(x, y)){
					x.degree++;
					y.degree++;
					x.possible=y;
					y.possible=x;
				}
			}
		}
		for(final Feature x:a){
			if(!x.exact && x.degree==1 && x.possible.degree==1){
				x.partner=x.possible;
				x.partner.partner=x;
			}
		}
		// Include unchanged anchors: retaining an old call while adding an overlapping
		// second call is still a split candidate, not a wholly unrelated gain.
		for(final Feature x:a){
			int reciprocal=0;
			final Group group=group(bc, x.seqid);
			for(int k=group.first(x.start); k<group.rows.size(); k++){
				final Feature y=group.rows.get(k);
				if(y.start>x.stop){break;}
				if(compatible(x, y) && overlap(x.start, x.stop, y.start, y.stop)>0){
					if(reciprocal(x, y)){reciprocal++;}
					if(x.cohort>=1 && x.cohort<=3 || !x.exact || !y.exact){
						reports.links.writer.print(reports.row.clear().append(binding.ordinal).tab().append(binding.id)
							.tab().append(x.line).tab().append(y.line).tab().append(overlap(x.start, x.stop, y.start, y.stop))
							.tab().append(reciprocal(x, y)).tab().append(x.exact).tab().append(y.exact).tab().append(x.partner==y).nl());
					}
				}
			}
			if(x.cohort>=1 && x.cohort<=3 && reciprocal>1){add(x, 17);}
		}
		for(final Feature y:b){
			if(y.cohort<1 || y.cohort>3){continue;}
			int reciprocal=0;
			final Group group=group(ac, y.seqid);
			for(int k=group.first(y.start); k<group.rows.size(); k++){
				final Feature x=group.rows.get(k);
				if(x.start>y.stop){break;}
				if(compatible(x, y) && reciprocal(x, y)){reciprocal++;}
			}
			if(reciprocal>1){add(y, 18);}
		}
		for(final Feature x:a){
			if(!x.exact){changed(binding, x, false, group(bc, x.seqid), ac, bc, reports);}
		}
		for(final Feature x:b){
			if(!x.exact){changed(binding, x, true, group(ac, x.seqid), ac, bc, reports);}
		}
		for(final Reference ref:binding.references){
			final int before=cover(ref, group(ac, ref.seqid));
			final int after=cover(ref, group(bc, ref.seqid));
			for(int pass=0; pass<2; pass++){
				final int scope=pass==0 ? 0 : ref.scored ? 1 : 2;
				final long[] c=referenceCounts[scope];
				c[0]++;
				if(before>0){c[1]++;}
				if(after>0){c[2]++;}
				if(before==0 && after>0){c[3]++;}
				if(before>0 && after==0){c[4]++;}
			}
			reports.coverage.writer.print(reports.row.clear().append(binding.ordinal).tab().append(binding.id).tab().append(ref.id)
				.tab().append(ref.seqid).tab().append(ref.start).tab().append(ref.stop).tab().append(ref.strand)
				.tab().append(ref.scored).tab().append(before).tab().append(after).nl());
		}
	}
	/** Complete within-path overlaps, not only changed loci. Each physical pair
	 * is emitted once even when both endpoints are clan calls. Multiplicity and
	 * opposite-strand overlaps stay explicit; these are not biological errors. */
	private void pathOverlaps(final Binding binding, final ArrayList<Feature> rows,
			final ObjectMap<String, Group> contigs, final int side, final Reports reports){
		final long[] counts=overlapCounts[side];
		for(final Feature anchor:rows){
			if(anchor.cohort<1 || anchor.cohort>3){continue;}
			final Group group=group(contigs, anchor.seqid);
			boolean affected=false;
			for(int i=group.first(anchor.start); i<group.rows.size(); i++){
				final Feature other=group.rows.get(i);
				if(other.start>anchor.stop){break;}
				if(other==anchor){continue;}
				final long bases=overlap(anchor.start, anchor.stop, other.start, other.stop);
				if(bases<1){continue;}
				affected=true;
				if(other.cohort>=1 && other.cohort<=3 && other.line<anchor.line){continue;}
				counts[0]++;
				counts[2]=Math.addExact(counts[2], bases);
				counts[3]=Math.max(counts[3], bases);
				reports.overlaps.writer.print(reports.row.clear().append(side==0 ? "BEFORE" : "AFTER")
					.tab().append(binding.ordinal).tab().append(binding.id).tab().append(anchor.seqid)
					.tab().append(anchor.line).tab().append(COHORTS[anchor.cohort]).tab().append(anchor.start).tab().append(anchor.stop).tab().append(anchor.strand)
					.tab().append(other.line).tab().append(COHORTS[other.cohort]).tab().append(other.type).tab().append(other.start).tab().append(other.stop).tab().append(other.strand)
					.tab().append(bases).nl());
			}
			if(affected){counts[1]++;}
		}
	}
	private void inventory(final ArrayList<Feature> rows, final ObjectMap<String, Group> groups, final int side){
		Feature previous=null;
		for(final Feature f:rows){
			add(f, side);
			if(previous!=null && previous.compareTo(f)==0){add(f, 13+side);}
			final Group group=groups.get(f.seqid);
			SpanCounts spans=group.spans.get(f.type);
			if(spans==null){
				spans=new SpanCounts();
				group.spans.put(f.type, spans);
			}
			// Two positive 31-bit coordinates and one strand bit fit in63 bits.
			final long span=((long)f.start<<32)|((long)f.stop<<1)|(f.strand=='-' ? 1L : 0L);
			if(!spans.add(f.cohort, span)){counts[f.cohort][15+side]++;}
			if(f.cohort>=1 && f.cohort<=3 && !spans.add(0, span)){counts[0][15+side]++;}
			previous=f;
		}
	}
	private void changed(final Binding binding, final Feature f, final boolean added, final Group opposite,
			final ObjectMap<String, Group> before, final ObjectMap<String, Group> after, final Reports reports){
		add(f, added ? 4 : 3);
		int overlaps=0;
		for(int k=opposite.first(f.start); k<opposite.rows.size(); k++){
			final Feature g=opposite.rows.get(k);
			if(g.start>f.stop){break;}
			if(compatible(f, g) && overlap(f.start, f.stop, g.start, g.stop)>0){overlaps++;}
		}
		if(overlaps==0){add(f, added ? 10 : 9);}
		if(f.partner!=null){
			final boolean coordinate=f.start!=f.partner.start || f.stop!=f.partner.stop;
			add(f, coordinate ? (added ? 6 : 5) : (added ? 8 : 7));
		}else if(overlaps>0){add(f, added ? 12 : 11);}
		int references=0;
		int scored=0;
		for(final Reference ref:binding.references){
			if(!ref.seqid.equals(f.seqid) || ref.strand!=f.strand || overlap(ref.start, ref.stop, f.start, f.stop)<1){continue;}
			references++;
			if(ref.scored){scored++;}
			reports.refLinks.writer.print(reports.row.clear().append(binding.ordinal).tab().append(binding.id).tab().append(added ? "ADDED" : "REMOVED")
				.tab().append(f.line).tab().append(ref.id).tab().append(ref.scored).tab().append(overlap(ref.start, ref.stop, f.start, f.stop)).nl());
		}
		final ByteBuilder row=reports.row.clear().append(binding.ordinal).tab().append(binding.id).tab().append(added ? "ADDED" : "REMOVED")
			.tab().append(f.line).tab().append(COHORTS[f.cohort]).tab();
		if(f.partner!=null){row.append(f.partner.line);}
		row.tab().append(f.degree).tab().append(overlaps).tab().append(references).tab().append(scored).tab().append(f.raw).nl();
		reports.changes.writer.print(row);
		if(f.cohort>=1 && f.cohort<=3){
			neighbours(binding, f, added, group(before, f.seqid).rows, "BEFORE", reports);
			neighbours(binding, f, added, group(after, f.seqid).rows, "AFTER", reports);
		}
	}

	/** All overlaps and nearest flanks on each strand, retaining ties. */
	private static void neighbours(final Binding binding, final Feature anchor, final boolean added,
			final ArrayList<Feature> rows, final String side, final Reports reports){
		final Feature[] left=new Feature[2];
		final Feature[] right=new Feature[2];
		for(final Feature f:rows){
			if(f==anchor){continue;}
			final int strand=f.strand=='+' ? 0 : 1;
			if(f.stop<anchor.start && (left[strand]==null || f.stop>left[strand].stop)){left[strand]=f;}
			if(f.start>anchor.stop && (right[strand]==null || f.start<right[strand].start)){right[strand]=f;}
		}
		for(final Feature f:rows){
			if(f==anchor){continue;}
			final int strand=f.strand=='+' ? 0 : 1;
			final long overlap=overlap(anchor.start, anchor.stop, f.start, f.stop);
			final String relation=f==anchor.partner ? "COUNTERPART" : overlap>0 ? "OVERLAP" : left[strand]!=null && f.stop==left[strand].stop ? "LEFT"
				: right[strand]!=null && f.start==right[strand].start ? "RIGHT" : null;
			if(relation!=null){reports.neighbours.writer.print(reports.row.clear().append(binding.ordinal).tab().append(binding.id)
				.tab().append(added ? "ADDED" : "REMOVED").tab().append(anchor.line).tab().append(side).tab().append(f.line)
				.tab().append(relation).tab().append(overlap).tab().append(f.exact).tab().append(f.raw).nl());}
		}
	}
	private static int cover(final Reference ref, final Group group){
		int count=0;
		for(int k=group.first(ref.start); k<group.rows.size(); k++){
			final Feature f=group.rows.get(k);
			if(f.start>ref.stop){break;}
			if(f.cohort>=1 && f.cohort<=3 && f.strand==ref.strand && overlap(f.start, f.stop, ref.start, ref.stop)>0){count++;}
		}
		return count;
	}
	private static ObjectMap<String, Group> contigs(final ArrayList<Feature> rows){
		final ObjectMap<String, Group> map=new ObjectMap<String, Group>(String.class, Group.class);
		final ArrayList<Group> groups=new ArrayList<Group>();
		for(final Feature f:rows){
			Group group=map.get(f.seqid);
			if(group==null){
				group=new Group();
				map.put(f.seqid, group);
				groups.add(group);
			}
			group.rows.add(f);
		}
		for(final Group group:groups){group.index();}
		return map;
	}
	private static Group group(final ObjectMap<String, Group> map, final String key){
		final Group result=map.get(key);
		return result==null ? EMPTY : result;
	}
	private void add(final Feature f, final int column){
		counts[f.cohort][column]++;
		if(f.cohort>=1 && f.cohort<=3){counts[0][column]++;}
	}
	private static boolean compatible(final Feature a, final Feature b){
		return a.seqid.equals(b.seqid) && a.strand==b.strand && a.type.equals(b.type)
			&& (a.cohort>=1 && a.cohort<=3 ? b.cohort>=1 && b.cohort<=3
				: a.cohort==b.cohort && a.kind.equals(b.kind));
	}
	private static boolean reciprocal(final Feature a, final Feature b){
		final long n=overlap(a.start, a.stop, b.start, b.stop);
		return n>0 && 2*n>=a.stop-a.start+1L && 2*n>=b.stop-b.start+1L;
	}
	private static long overlap(final int a, final int b, final int c, final int d){return Math.max(0L, Math.min(b, d)-(long)Math.max(a, c)+1);}
	private static int integer(final LineParser1 parser, final int column){
		parser.setBounds(column);
		return NcrnaNumbers.unsignedInt(parser.line(), parser.a(), parser.b());
	}
	/** Validate the retained GFF score without native integer-accumulator overflow.
	 * Ordinary short decimal fields stay on borrowed bytes; long or exponent
	 * fields use the JDK parser. No numeric score is used for matching. */
	private static void score(final LineParser1 parser){
		if(parser.termEquals(".", 5)){return;}
		parser.setBounds(5);
		final byte[] line=parser.line();
		final int end=parser.b();
		int i=parser.a();
		final boolean plus=i<end && line[i]=='+';
		if(i<end && (line[i]=='-' || line[i]=='+')){i++;}
		int digits=0;
		boolean dot=false;
		while(i<end){
			final byte c=line[i];
			if(c>='0' && c<='9'){
				digits++;
				i++;
			}else if(c=='.' && !dot){
				dot=true;
				i++;
			}
			else{break;}
		}
		require(digits>0, "GFF score must be a decimal number or an explicit dot");
		final boolean exponent=i<end && (line[i]=='e' || line[i]=='E');
		if(exponent){
			i++;
			if(i<end && (line[i]=='-' || line[i]=='+')){i++;}
			final int from=i;
			while(i<end && line[i]>='0' && line[i]<='9'){i++;}
			require(i>from, "GFF score exponent requires decimal digits");
		}
		require(i==end, "GFF score contains invalid decimal syntax");
		final float value=plus || exponent || digits>9 ? Float.parseFloat(parser.parseString(5)) : parser.parseFloat(5);
		require(Float.isFinite(value), "Computed caller score must be finite");
	}
	private Output output(final String name){return new Output(new File(out, name).getPath());}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Cursor implements AutoCloseable{
		Cursor(final String path){
			this.path=path;
			input=new Input(path);
			try{
				pending=next();
				while(pending!=null && !boundary(pending)){
					require(pending.length>0 && pending[0]=='#', "Feature precedes its input ordinal: "+path);
					pending=next();
				}
			}catch(RuntimeException | Error failure){
				try{input.close();}catch(RuntimeException close){failure.addSuppressed(close);}
				throw failure;
			}
		}
		ArrayList<Feature> input(final int ordinal){
			require(pending!=null && boundary(pending), "Missing GFF input ordinal "+ordinal);
			parser.set(pending);
			require(parser.terms()==2 && integer(parser, 1)==ordinal, "Input ordinals must be contiguous from zero");
			final ArrayList<Feature> result=new ArrayList<Feature>();
			pending=next();
			while(pending!=null && !boundary(pending)){
				if(pending.length>0 && pending[0]!='#'){
					try{
						parser.set(pending);
						require(parser.terms()==9, "Saved caller feature needs nine fields");
						if(seqid==null || !parser.termEquals(seqid, 0)){seqid=parser.parseString(0);}
						if(type==null || !parser.termEquals(type, 2)){
							type=parser.termEquals("CDS", 2) ? "CDS" : parser.termEquals("RNA", 2) ? "RNA"
								: parser.termEquals("tRNA", 2) ? "tRNA" : parser.parseString(2);
						}
						result.add(new Feature(parser, line, seqid, type));
					}
					catch(IllegalArgumentException failure){throw new IllegalArgumentException("Invalid saved GFF at "+path+":"+line+": "+failure.getMessage(), failure);}
				}
				pending=next();
			}
			return result;
		}
		byte[] next(){
			final byte[] row=input.reader.nextLine();
			if(row!=null){
				line++;
				require(!shared.Tools.startsWith(row, "##FASTA"), "Caller GFF cannot contain sequence payload");
			}
			return row;
		}
		static boolean boundary(final byte[] row){return shared.Tools.startsWith(row, "##input_ordinal0\t");}
		@Override
		public void close(){input.close();}
		final Input input;
		final String path;
		final LineParser1 parser=new LineParser1('\t');
		byte[] pending;
		String seqid;
		String type;
		long line;
	}
	private static final class Feature implements Comparable<Feature>{
		Feature(final LineParser1 p, final long line_, final String seqid_, final String type_){
			require(p.terms()==9 && p.length(0)>0 && p.length(2)>0 && (p.termEquals("+", 6) || p.termEquals("-", 6)), "Saved caller feature needs nine fields and an explicit strand");
			start=integer(p, 3);
			stop=integer(p, 4);
			require(start>0 && stop>=start, "Saved caller coordinates must be positive closed intervals");
			score(p);
			seqid=seqid_;
			type=type_;
			strand=p.termEquals("+", 6) ? '+' : '-';
			int member=0;
			p.setBounds(8);
			final byte[] bytes=p.line();
			final int from=p.a();
			final int to=p.b();
			final boolean clan=shared.Tools.indexOf(bytes, "model:cobalamin_", from)>=0;
			if(clan){
				for(int i=1; i<=3; i++){
					if(shared.Tools.indexOf(bytes, MEMBER_PREFIXES[i-1], from)>=0){
						require(member==0, "One feature cannot belong to two cobalamin members");
						member=i;
					}
				}
			}
			require(member>0 || !clan, "Unknown cobalamin member must not silently become OTHER");
			require(member==0 || type.equals("RNA"), "A cobalamin model must belong to an RNA feature");
			cohort=member>0 ? member : type.equals("CDS") ? 4 : 5;
			if(cohort==5){
				final int comma=shared.Tools.indexOf(bytes, (byte)',', from);
				kind=type+'\t'+new String(bytes, from, (comma<0 ? to : comma)-from, StandardCharsets.UTF_8);
			}else{kind=type;}
			// ByteFile returns an owned line array; its parser only borrows it.
			// Retain that array for this input instead of copying every GFF row.
			raw=p.line();
			line=line_;
		}
		@Override
		public int compareTo(final Feature other){return GffCallSetDelta.compare(raw, other.raw);}
		final byte[] raw;
		final String seqid;
		final String type;
		final String kind;
		final int start;
		final int stop;
		final int cohort;
		final char strand;
		final long line;
		boolean exact;
		int degree;
		Feature possible;
		Feature partner;
	}
	private static final class Reference{
		Reference(final LineParser1 p){
			id=p.parseString(1);
			seqid=p.parseString(2);
			start=integer(p, 3);
			stop=integer(p, 4);
			require(!id.isEmpty() && !seqid.isEmpty() && start>0 && stop>=start
				&& (p.termEquals("+", 5) || p.termEquals("-", 5)) && (p.termEquals("true", 6) || p.termEquals("false", 6)), "Frozen reference geometry, strand and scored status must be explicit");
			strand=p.termEquals("+", 5) ? '+' : '-';
			scored=p.termEquals("true", 6);
		}
		final String id;
		final String seqid;
		final int start;
		final int stop;
		final char strand;
		final boolean scored;
	}
	private static final class Binding{
		Binding(final int ordinal_, final String id_, final String fasta_){
			ordinal=ordinal_;
			id=id_;
			fasta=fasta_;
		}
		final int ordinal;
		final String id;
		final String fasta;
		final ArrayList<Reference> references=new ArrayList<Reference>();
	}
	private static final class Group{
		/** Sort once per contig. Prefix maxima retain long enclosing intervals;
		 * searching only starts would silently miss those overlaps. */
		void index(){
			Collections.sort(rows, POSITION);
			maxStop=new int[rows.size()];
			int maximum=0;
			for(int i=0; i<rows.size(); i++){
				maximum=Math.max(maximum, rows.get(i).stop);
				maxStop[i]=maximum;
			}
		}
		int first(final int start){
			assert(maxStop.length==rows.size()) : "Every contig is indexed before overlap queries; missing suffixes would silently drop overlaps";
			int low=0;
			int high=maxStop.length;
			while(low<high){
				final int middle=(low+high)>>>1;
				if(maxStop[middle]<start){low=middle+1;}else{high=middle;}
			}
			return low;
		}
		final ArrayList<Feature> rows=new ArrayList<Feature>();
		final ObjectMap<String, SpanCounts> spans=new ObjectMap<String, SpanCounts>(String.class, SpanCounts.class);
		int[] maxStop=new int[0];
	}
	private static final class SpanCounts{
		boolean add(final int cohort, final long span){
			assert(cohort>=0 && cohort<seen.length && span>=0) : "Only validated positive coordinates and declared cohorts enter duplicate accounting";
			LongHashSet set=seen[cohort];
			if(set==null){
				set=new LongHashSet(16);
				seen[cohort]=set;
			}
			return set.add(span);
		}
		final LongHashSet[] seen=new LongHashSet[6];
	}
	private static final class Reports{
		Reports(final Output a, final Output b, final Output c, final Output d, final Output e, final Output f){
			changes=a;
			links=b;
			neighbours=c;
			coverage=d;
			refLinks=e;
			overlaps=f;
		}
		void headers(){
			changes.writer.println("input_ordinal0\tinput_id\tchange\tsource_line\tcohort\tpartner_line_if_unique\treciprocal_candidates\tall_counterpart_overlaps\treference_overlaps\tscored_reference_overlaps\t"+GFF);
			links.writer.println("input_ordinal0\tinput_id\tbefore_line\tafter_line\toverlap_bp\treciprocal_half\tbefore_exact\tafter_exact\tunique_unmatched_pair");
			neighbours.writer.println("input_ordinal0\tinput_id\tanchor_change\tanchor_line\tneighbour_side\tneighbour_line\trelation\toverlap_bp\texact_in_both\t"+GFF);
			coverage.writer.println("input_ordinal0\tinput_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\tclan_overlaps_before\tclan_overlaps_after");
			refLinks.writer.println("input_ordinal0\tinput_id\tchange\tsource_line\treference_id\tscored\toverlap_bp");
			overlaps.writer.println("side\tinput_ordinal0\tinput_id\tseqid\tclan_line\tclan_member\tclan_start\tclan_stop\tclan_strand\tneighbour_line\tneighbour_cohort\tneighbour_type\tneighbour_start\tneighbour_stop\tneighbour_strand\toverlap_bp");
		}
		final Output changes;
		final Output links;
		final Output neighbours;
		final Output coverage;
		final Output refLinks;
		final Output overlaps;
		final ByteBuilder row=new ByteBuilder();
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String before;
	private String after;
	private String beforeConfig;
	private String afterConfig;
	private String sources;
	private String genomes;
	private String reference;
	private String out;
	private boolean legacyNames;
	private int expected=-1;
	private final ArrayList<Binding> order=new ArrayList<Binding>();
	private final ObjectMap<String, Binding> byId=new ObjectMap<String, Binding>(String.class, Binding.class);
	private final long[][] counts=new long[6][19];
	private final long[][] referenceCounts=new long[3][5];
	private final long[][] overlapCounts=new long[2][4];
	private static final Group EMPTY=new Group();
	private static final Comparator<Feature> POSITION=(a, b)->{
		final int start=Integer.compare(a.start, b.start);
		return start!=0 ? start : Integer.compare(a.stop, b.stop);
	};
	private static final String[] COHORTS={"CLAN", "RF00174", "RF01482", "RF01689", "CDS", "OTHER"};
	private static final String[] MEMBER_PREFIXES={"model:cobalamin_RF00174_", "model:cobalamin_RF01482_", "model:cobalamin_RF01689_"};
	private static final String GFF="seqid\tsource\ttype\tstart\tstop\tscore\tstrand\tphase\tattributes";
}
