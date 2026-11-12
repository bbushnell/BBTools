package prok;

import java.util.Arrays;
import java.util.Comparator;

import parse.LineParser1;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaReferenceGeometry.InputState;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;
import structures.DoubleList;

/** Describes annotation-concordant observations from previously validated row links.
 * A source link must identify the exact input/native table/line and explicitly
	 * declare same-strand overlap and a matching accession or biological family.
	 * Literal tRNA/tmRNA feature types with explicitly absent RFAM claims provide
	 * a separately labeled type-only family basis; they never identify a subtype.
 * RNA type, RFAM agreement and warnings are read from that same source-link row;
 * event-wide OR flags are never used as a substitute. The upstream strata producer
 * established literal-contig and positive inclusive-overlap geometry. Its saved
 * source links and original evidence must be supplied as one accepted generation.
 *
 * Physical intervals are keyed by input, contig, strand, bounds and model. Repeated
 * annotation links or identical native observations do not weight the distribution.
 * Conflicting native observations at identical bounds retain one anchor but no
 * chosen score. Source warnings are the union of matching source-row claims, not
 * inferred RNA properties. Source-kind counts are nonexclusive. Exact-accession
 * support takes precedence over family-only support, but neither establishes
 * independent truth, a subtype decision, an FPR or an admission threshold.
 * All mutable state is owned by one calling thread; numeric data use native lists.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceSourceAnchors{

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	static void process(final String input, final String strata, final String output,
		final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceSourceAnchors report=new NcrnaReferenceSourceAnchors(input, catalog, roster);
		report.read(strata);
		report.write(output);
	}
	private NcrnaReferenceSourceAnchors(final String input, final NcrnaEvidenceCatalog catalog_, final Roster roster_){
		assert(TYPE!=0 && (TYPE&(WARNING_MASK|63))==0) : "Type-only support must not overlap RFAM, source-kind or warning bits";
		catalog=catalog_;
		roster=roster_;
		geometry=NcrnaReferenceGeometry.load(input, catalog, roster);
		flags=new int[geometry.hits().length];
		links=new int[flags.length];
		firstExact=new long[flags.length];
		firstFamily=new long[flags.length];
		firstType=new long[flags.length];
		stats=new Stats[roster.queried()][BASES.length];
		for(final Stats[] row:stats){
			for(int i=0; i<row.length; i++){row[i]=new Stats(catalog.size());}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Qualified links        ----------------*/
	/*--------------------------------------------------------------*/

	private void read(final String prefix){
		linkFile=prefix+".source_links.tsv.gz";
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(prefix+".strata_summary.tsv")){
			RrnaResourceIO.header(in.reader, STRATA_HEADER, prefix);
			final byte[] row=in.reader.nextLine();
			require(row!=null, "Saved source strata need a completed summary");
			p.set(row);
			require(p.terms()==8 && p.parseInt(0)==geometry.profiles().observed()
				&& p.parseInt(1)==catalog.size() && p.parseLong(2)==flags.length
				&& p.parseInt(3)==geometry.groups() && p.parseLong(4)==geometry.profiles().retained()
				&& p.parseLong(5)==geometry.profiles().historicalPairs()
				&& p.termEquals(geometry.profiles().evidenceState(), 6)
				&& p.termEquals("STRATA_COMPLETE_NOT_ADJUDICATED", 7) && in.reader.nextLine()==null,
				"Source-link population must match the supplied saved geometry");
		}
		InputState input=null;
		InputState gffInput=null;
		byte[] gffPath=null;
		final ByteBuilder encoded=new ByteBuilder();
		try(Input in=new Input(linkFile)){
			RrnaResourceIO.header(in.reader, LINK_HEADER, linkFile);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				require(p.terms()==10 && p.parseLong(2)>0, "Saved source links need ten fields and a physical source line");
				if(input==null || !p.termEquals(input.genome.id, 0)){input=geometry.input(p.parseString(0));}
				require(input!=null && input.selected,
					"Source link must retain its exact observed input and native table");
				if(gffInput!=input || !same(p, 1, gffPath)){
					boolean claimed=false;
					for(int m=0; m<input.genome.memberships(); m++){
						BuildNcrnaReference.field(encoded.clear(), input.genome.membership(m).gff);
						p.setBounds(1);
						if(p.b()-p.a()!=encoded.length()){continue;}
						boolean equal=true;
						for(int j=0; equal && j<encoded.length(); j++){equal=p.line()[p.a()+j]==encoded.array[j];}
						claimed|=equal;
					}
					require(claimed, "Source GFF must be declared by this exact catalog input");
					gffInput=input;
					gffPath=p.parseByteArray(1);
				}
				final int index=input.find(p, 3, p.parseLong(4));
				require(index>=0, "Source link names no saved native event");
				final Hit hit=geometry.hits()[index];
				int kind=-1;
				for(int i=0; i<KINDS.length; i++){if(p.termEquals(KINDS[i], 6)){kind=i;}}
				require(kind>=0, "Unknown source-kind qualification");
				linkRows++;
				if(!p.termEquals("SAME_STRAND_OVERLAP", 7)){
					require(p.termEquals("OPPOSITE_STRAND_OVERLAP", 7) || p.termEquals("UNASSESSABLE", 7),
						"Source links require an explicit spatial association");
					continue;
				}
				require(hit.valid(), "Unassessable native coordinates cannot supply annotation concordance");
				final boolean exact=token(p, 8, "SOURCE_SAME_ACCESSION");
				final boolean family=token(p, 8, "SOURCE_SAME_FAMILY_OTHER_ACCESSION");
				final String hitFamily=roster.query(hit.query).scope.family;
				final boolean typed=kind==0 && p.termEquals("NO_OVERLAPPING_SOURCE_RFAM", 8)
					&& (p.termEquals("tRNA", 5) && hitFamily.equals("tRNA") || p.termEquals("tmRNA", 5) && hitFamily.equals("tmRNA"));
				if(!exact && !family && !typed){continue;}
				int word=(exact ? EXACT : 0)|(family ? FAMILY : 0)|(typed ? TYPE : 0)|(1<<(kind+2));
				for(int i=0; i<WARNINGS.length; i++){
					if(token(p, 8, WARNINGS[i]) || token(p, 9, WARNINGS[i])){word|=1<<(i+6);}
				}
				flags[index]|=word;
				links[index]=Math.addExact(links[index], 1);
				qualifiedLinks++;
				if(exact && firstExact[index]==0){firstExact[index]=line;}
				if(family && firstFamily[index]==0){firstFamily[index]=line;}
				if(typed && firstType[index]==0){firstType[index]=line;}
			}
		}
		require(linkRows==geometry.profiles().historicalPairs(), "Source-link rows must conserve the accepted population");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Physical intervals     ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		final Hit[] ordered=geometry.hits().clone();
		Arrays.sort(ordered, ORDER);
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(prefix+SUFFIXES[0])){
			out.writer.println("input_id\ttarget\tstrand\tstart1\tstop1\taccession\tbasis\tnative_events\tsupport_links\tlink_file\tfirst_exact_link_line\tfirst_family_link_line\tfirst_type_link_line\tsource_kinds\tsource_warnings\tscore_observations\tbits\tGA_margin_bits\tmodel_coverage\tnative_truncation\tleft_edge\tright_edge\tstatus");
			for(int i=0; i<ordered.length; ){
				final Hit first=ordered[i];
				int end=i+1, word=flags[first.index];
				long count=links[first.index], exact=firstExact[first.index], family=firstFamily[first.index];
				long typed=firstType[first.index];
				boolean consistent=true;
				while(end<ordered.length && sameInterval(first, ordered[end])){
					final Hit hit=ordered[end++];
					word|=flags[hit.index];
					count+=links[hit.index];
					exact=firstNonzero(exact, firstExact[hit.index]);
					family=firstNonzero(family, firstFamily[hit.index]);
					typed=firstNonzero(typed, firstType[hit.index]);
					consistent&=sameObservation(first, hit);
				}
				physicalIntervals++;
				if((word&(EXACT|FAMILY|TYPE))!=0){
					final int basis=((word&EXACT)!=0 ? 0 : (word&FAMILY)!=0 ? 2 : 4)+((word&WARNING_MASK)==0 ? 0 : 1);
					final Stats s=stats[first.query][basis];
					final Model m=roster.query(first.query);
					s.anchors++;
					s.events+=end-i;
					s.supportLinks+=count;
					supportedIntervals++;
					if(!s.inputs[first.input.index]){
						s.inputs[first.input.index]=true;
						s.genomes++;
					}
					for(int k=0; k<KINDS.length; k++){if((word&(1<<(k+2)))!=0){s.kinds[k]++;}}
					if(end-i>1){s.duplicates++;}
					if(!consistent){s.conflicts++;}
					else{
						assert(first.valid()) : "Concordant support requires a valid native interval, and duplicate observations agree";
						s.margins.add(first.bits-m.ga);
						s.coverage.add((first.to-first.from+1)/(double)m.length);
						if(first.full()){
							s.full++;
							s.fullMargins.add(first.bits-m.ga);
						}
						if(first.left || first.right){s.edges++;}
					}
					b.clear().append(first.input.genome.id).tab().append(first.target.bytes).tab().append(first.strand).tab();
					b.append(first.start).tab().append(first.stop).tab().append(m.accession).tab().append(BASES[basis]).tab();
					b.append(end-i).tab().append(count).tab();
					BuildNcrnaReference.field(b, linkFile).tab();
					witness(b, exact).tab();
					witness(b, family).tab();
					witness(b, typed).tab();
					kindNames(b, word).tab();
					warningNames(b, word).tab().append(consistent ? 1 : 0).tab();
					if(consistent){
						b.append(first.bits, 6).tab().append(first.bits-m.ga, 6).tab().append((first.to-first.from+1)/(double)m.length, 6).tab();
						b.append(first.truncation).tab().append(first.left).tab().append(first.right);
					}else{b.append("NA\tNA\tNA\tNA\tNA\tNA");}
					out.writer.print(b.append("\tANNOTATION_CONCORDANT_NOT_REVIEWED_TRUTH\n"));
				}
				i=end;
			}
		}
		writeModels(prefix+SUFFIXES[1]);
		try(Output out=new Output(prefix+SUFFIXES[2])){
			out.writer.println("inputs_observed\tcatalog_inputs\tnative_events\tphysical_model_intervals\tsource_link_rows\tqualified_links\tconcordant_intervals\tevidence_status\tstatus");
			b.clear().append(geometry.profiles().observed()).tab().append(catalog.size()).tab().append(ordered.length).tab();
			b.append(physicalIntervals).tab().append(linkRows).tab().append(qualifiedLinks).tab().append(supportedIntervals).tab();
			out.writer.print(b.append(geometry.profiles().evidenceState()).append("\tANNOTATION_CONCORDANCE_COMPLETE_NOT_ADMISSION\n"));
		}
	}
	private void writeModels(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			out.writer.println("accession\trelease\tmodel_sha80\tGA_bits\tmodel_length\tbasis\tphysical_anchors\tdistinct_inputs\tnative_events\tsupport_links\trepeated_native_intervals\tconflicting_native_intervals\tscore_observations\tfull_untruncated\tedge_observations\trna_supported\tcomponent_supported\tcontext_supported\tother_supported\tmargin_min\tmargin_p05\tmargin_median\tmargin_p95\tmargin_max\tcoverage_min\tcoverage_p05\tcoverage_median\tcoverage_p95\tcoverage_max\tfull_margin_min\tfull_margin_p05\tfull_margin_median\tfull_margin_p95\tfull_margin_max\tstatus");
			long total=0, linkTotal=0;
			for(int q=0; q<stats.length; q++){
				final Model m=roster.query(q);
				for(int basis=0; basis<BASES.length; basis++){
					final Stats s=stats[q][basis];
					require(s.anchors==s.conflicts+s.margins.size, "Each physical anchor needs one score observation or an explicit conflict");
					total+=s.anchors;
					linkTotal+=s.supportLinks;
					b.clear().append(m.accession).tab().append(m.release).tab().append(m.sha80).tab().append(m.ga, 6).tab().append(m.length).tab();
					b.append(BASES[basis]).tab().append(s.anchors).tab().append(s.genomes).tab().append(s.events).tab().append(s.supportLinks).tab();
					b.append(s.duplicates).tab().append(s.conflicts).tab().append(s.margins.size).tab().append(s.full).tab().append(s.edges);
					for(final long value:s.kinds){b.tab().append(value);}
					NcrnaReferenceAnchors.quantiles(b, s.margins);
					NcrnaReferenceAnchors.quantiles(b, s.coverage);
					NcrnaReferenceAnchors.quantiles(b, s.fullMargins);
					out.writer.print(b.append("\tOBSERVED_SOURCE_CONCORDANCE_NOT_FPR\n"));
				}
			}
			require(total==supportedIntervals && linkTotal==qualifiedLinks, "Physical grouping must conserve all qualifying source links without weighting loci by link count");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Byte helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean token(final LineParser1 p, final int field, final String name){
		p.setBounds(field);
		final byte[] row=p.line();
		for(int start=p.a(); start<p.b(); ){
			int end=start;
			while(end<p.b() && row[end]!=';'){end++;}
			boolean equal=end-start==name.length();
			for(int i=0; equal && i<name.length(); i++){equal=row[start+i]==name.charAt(i);}
			if(equal){return true;}
			start=end+1;
		}
		return false;
	}
	private static boolean same(final LineParser1 p, final int field, final byte[] bytes){
		p.setBounds(field);
		if(bytes==null || p.b()-p.a()!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(bytes[i]!=p.line()[p.a()+i]){return false;}}
		return true;
	}
	private static boolean sameInterval(final Hit a, final Hit b){
		return a.input==b.input && a.target.equals(b.target) && a.strand==b.strand
			&& a.start==b.start && a.stop==b.stop && a.query==b.query;
	}
	private static boolean sameObservation(final Hit a, final Hit b){
		return a.bits==b.bits && a.from==b.from && a.to==b.to && a.truncation==b.truncation
			&& a.range==b.range && a.left==b.left && a.right==b.right && a.ga==b.ga;
	}
	private static long firstNonzero(final long a, final long b){return a==0 ? b : b==0 ? a : Math.min(a, b);}
	private static ByteBuilder witness(final ByteBuilder b, final long line){return line==0 ? b.append("NA") : b.append(line);}
	private static ByteBuilder kindNames(final ByteBuilder b, final int word){
		assert((word&60)!=0) : "Every concordant interval retains at least one directly matching source-row kind";
		boolean first=true;
		for(int i=0; i<KINDS.length; i++){
			if((word&(1<<(i+2)))==0){continue;}
			if(!first){b.append(';');}
			b.append(KINDS[i]);
			first=false;
		}
		return b;
	}
	private static ByteBuilder warningNames(final ByteBuilder b, final int word){
		if((word&WARNING_MASK)==0){return b.append("NONE");}
		boolean first=true;
		for(int i=0; i<WARNINGS.length; i++){
			if((word&(1<<(i+6)))==0){continue;}
			if(!first){b.append(';');}
			b.append(WARNINGS[i]);
			first=false;
		}
		return b;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             State            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Stats{
		Stats(final int count){inputs=new boolean[count];}
		final boolean[] inputs;
		final long[] kinds=new long[4];
		final DoubleList margins=new DoubleList(1), coverage=new DoubleList(1), fullMargins=new DoubleList(1);
		long anchors, genomes, events, supportLinks, duplicates, conflicts, full, edges;
	}
	private final NcrnaEvidenceCatalog catalog;
	private final Roster roster;
	private final NcrnaReferenceGeometry geometry;
	private final int[] flags, links;
	private final long[] firstExact, firstFamily, firstType;
	private final Stats[][] stats;
	private String linkFile;
	private long linkRows, qualifiedLinks, physicalIntervals, supportedIntervals;
	static final String[] SUFFIXES={".source_anchor_loci.tsv.gz", ".source_anchor_models.tsv", ".source_anchor_summary.tsv"};
	private static final int EXACT=1, FAMILY=2;
	private static final String[] BASES={"EXACT_ACCESSION_CLEAR", "EXACT_ACCESSION_FLAGGED", "FAMILY_ONLY_CLEAR", "FAMILY_ONLY_FLAGGED", "TYPE_FAMILY_ONLY_CLEAR", "TYPE_FAMILY_ONLY_FLAGGED"};
	private static final String[] KINDS={"RNA_ROW", "COMPONENT_ROW", "GENE_OR_TRANSCRIPT_CONTEXT", "OTHER_ROW"};
	private static final String[] WARNINGS={"SOURCE_PSEUDO_TRUE", "SOURCE_PSEUDOGENE_ATTRIBUTE", "SOURCE_PSEUDOGENE_BIOTYPE", "SOURCE_PARTIAL_TRUE",
		"SOURCE_EXCEPTION_ATTRIBUTE", "SOURCE_UNRECOGNIZED_BOOLEAN", "SOURCE_INCOMPATIBLE_ID", "SOURCE_DANGLING_PARENT", "SOURCE_UNKNOWN_RFAM",
		"SOURCE_OTHER_KNOWN_FAMILY", "SOURCE_RELATED_CLASS", "UNASSESSABLE_RANGE_OR_SOURCE", "STRUCTURAL_SOURCE_CONTEXT"};
	private static final int WARNING_MASK=((1<<WARNINGS.length)-1)<<6;
	private static final int TYPE=1<<(6+WARNINGS.length);
	private static final String LINK_HEADER="input_id\tgff\tline\tnative_table\tnative_line\tsource_type\tsource_kind\tassociation\trfam_comparison\tsource_flags";
	private static final String STRATA_HEADER="inputs_observed\tcatalog_inputs\traw_hit_events\treview_groups\tretained_source_rows\thistorical_model_pairs\tevidence_status\tstatus";
	private static final Comparator<Hit> ORDER=(a, b)->{
		int c=Integer.compare(a.input.index, b.input.index);
		if(c==0){c=Tools.compare(a.target.bytes, b.target.bytes);}
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		if(c==0){c=Integer.compare(a.query, b.query);}
		if(c==0){c=Tools.compare(a.table, b.table);}
		return c==0 ? Long.compare(a.line, b.line) : c;
	};
}
