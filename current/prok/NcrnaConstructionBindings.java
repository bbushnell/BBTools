package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Membership;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Binds the original FMN construction population to the two discovery cohorts.
 * Reads metadata only: neither a path match nor an assembly claim proves sequence
 * equivalence. Keeps every original seed input, including unassigned members,
 * and preserves every matching cohort/input rather than selecting an arbitrary
 * alias. This is input binding for an eligibility ledger, not admission or a split.
 * One owner constructs and writes this object; no worker mutates its maps.
 * @author Ganyu
 */
public final class NcrnaConstructionBindings{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionBindings.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String genomes=null, loci=null, assignments=null, evalCensus=null, evalMembers=null, newCensus=null, newMembers=null;
			int expectedGenomes=-1, expectedMembers=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("genomes")){genomes=b;}
				else if(a.equals("loci")){loci=b;}
				else if(a.equals("assignments")){assignments=b;}
				else if(a.equals("evalcensus")){evalCensus=b;}
				else if(a.equals("evalmembers")){evalMembers=b;}
				else if(a.equals("newcensus")){newCensus=b;}
				else if(a.equals("newmembers")){newMembers=b;}
				else if(a.equals("expectedgenomes")){expectedGenomes=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmembers")){expectedMembers=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown construction binding parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions, "Construction identity checks require native assertions enabled");
			require(genomes!=null && loci!=null && assignments!=null && evalCensus!=null && evalMembers!=null
				&& newCensus!=null && newMembers!=null && parser.out1!=null && expectedGenomes>0 && expectedMembers>0,
				"Supply both catalogs, the three original tables, expected counts and out=fresh_prefix");
			require(parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Complete original populations cannot be sampled or overwritten");
			final String[] files={genomes, loci, assignments, evalCensus, evalMembers, newCensus, newMembers};
			require(Tools.testInputFiles(false, true, files), "All declared metadata inputs must be readable");
			require(Tools.testForDuplicateFiles(true, files), "Original tables and cohort inputs must have distinct paths");
			final NcrnaConstructionBindings result=load(genomes, loci, assignments,
				NcrnaEvidenceCatalog.load(evalCensus, evalMembers), NcrnaEvidenceCatalog.load(newCensus, newMembers),
				expectedGenomes, expectedMembers);
			result.write(parser.out1);
			pp.outstream.println("NCRNA_CONSTRUCTION_BINDINGS_COMPLETE_NOT_ELIGIBILITY members="+result.members()+" kept="+result.kept());
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Joins complete original tables and existing path-based catalogs.
	 * No FASTA, GFF or covariance-model payload is opened here. */
	public static NcrnaConstructionBindings load(final String genomes, final String loci,
			final String assignments, final NcrnaEvidenceCatalog evaluation,
			final NcrnaEvidenceCatalog discovery, final int expectedGenomes, final int expectedMembers){
		require(evaluation!=null && discovery!=null && expectedGenomes>0 && expectedMembers>0,
			"Complete evaluation/discovery catalogs and explicit original population counts are required");
		final NcrnaConstructionBindings result=new NcrnaConstructionBindings();
		result.catalog(evaluation, true);
		result.catalog(discovery, false);
		result.genomes(genomes, expectedGenomes);
		result.loci(loci, expectedMembers);
		result.assignments(assignments);
		return result;
	}

	private NcrnaConstructionBindings(){}

	/*--------------------------------------------------------------*/
	/*----------------       Cohort Bindings        ----------------*/
	/*--------------------------------------------------------------*/

	private void catalog(final NcrnaEvidenceCatalog catalog, final boolean evaluation){
		assert(catalog.size()>0) : "A binding cohort must retain the complete nonempty catalog denominator";
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			final Binding binding=new Binding(evaluation, genome.id);
			path(genome.fasta, binding);
			for(int j=0; j<genome.memberships(); j++){
				final Membership member=genome.membership(j);
				path(member.suppliedFasta, binding);
				if(evaluation && !member.assembly.equals("NA")){evaluationAssemblies.add(member.assembly);}
				if(!member.assembly.equals("NA") && evaluationAssemblies.contains(member.assembly)){binding.evaluationAssembly=true;}
			}
		}
	}

	private void path(final String path, final Binding binding){
		require(path.startsWith("/"), "Catalog path claims must remain absolute");
		PathBindings values=paths.get(path);
		if(values==null){
			values=new PathBindings();
			paths.put(path, values);
		}
		for(final Binding old:values.bindings){
			if(old.evaluation==binding.evaluation && old.input.equals(binding.input)){return;}
		}
		values.bindings.add(binding);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Original Members       ----------------*/
	/*--------------------------------------------------------------*/

	private void genomes(final String path, final int expected){
		long members=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GENOMES_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11, "Original construction roster requires eleven fields");
				final int index=p.parseInt(0), count=p.parseInt(8);
				final String name=p.parseString(1), fasta=p.parseString(7);
				require(index==originalGenomes.size() && count>=0 && !name.isEmpty() && !byGenome.contains(name),
					"Original roster order, unique genome names and member counts must be preserved");
				final PathBindings matches=paths.get(fasta);
				require(matches!=null && !matches.bindings.isEmpty(), "Original genome has no literal cohort path binding: "+fasta);
				final OriginalGenome genome=new OriginalGenome(name, p.parseString(2), p.parseString(3),
					p.parseString(4), p.parseString(5), p.parseString(6), fasta, count, matches);
				require(genome.gff.startsWith("/") && !genome.assembly.isEmpty() && !genome.taxid.isEmpty(),
					"Original source identity must retain its paired GFF, assembly and taxonomy claims");
				genome.evaluationAssembly=!genome.assembly.equals("NA") && evaluationAssemblies.contains(genome.assembly);
				for(final Binding binding:matches.bindings){genome.evaluationAssembly|=binding.evaluationAssembly;}
				byGenome.put(name, genome);
				originalGenomes.add(genome);
				members=Math.addExact(members, count);
			}
		}
		require(originalGenomes.size()==expected, "Original genome roster differs from the declared complete population");
		expectedFromRoster=members;
	}

	private void loci(final String path, final int expected){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, LOCI_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==21, "Original locus ledger requires all twenty-one provenance fields");
				final String id=p.parseString(0), name=p.parseString(1);
				final OriginalGenome genome=byGenome.get(name);
				require(genome!=null && id.startsWith("RF00050_") && !byLocus.contains(id),
					"Each original FMN locus requires one unique identity and an original roster owner: "+id);
				require(p.termEquals(genome.assembly, 2) && p.termEquals(genome.taxid, 3)
					&& p.termEquals(genome.domain, 4) && p.termEquals(genome.phylum, 5),
					"Locus identity claims must match their original genome roster: "+id);
				final int length=p.parseInt(9), contigLength=p.parseInt(15);
				require(!p.termEquals("", 6) && (p.termEquals("+", 7) || p.termEquals("-", 7)) && length>0 && contigLength>=length,
					"Original loci need explicit contig, strand and positive bounded lengths: "+id);
				p.setBounds(8);
				final int intervals=checkIntervals(p.line(), p.a(), p.b(), length, contigLength);
				require(p.parseLong(10)>0 && !p.termEquals("", 11) && !p.termEquals("", 20),
					"Original source lines and feature IDs must survive the construction join: "+id);
				bool(12);
				bool(13);
				bool(16);
				final int undefined=p.parseInt(17);
				require(undefined>=0 && undefined<=length && p.parseInt(18)>=0 && p.parseInt(19)>=0,
					"Original ambiguity and clipped-flank counts must be nonnegative and bounded: "+id);
				final String status=p.parseString(14);
				final String expectedStatus=p.termEquals("true", 13) ? "PSEUDOGENE" : p.termEquals("true", 12)
					? "SOURCE_PARTIAL" : intervals>1 ? "MULTIPART_UNRESOLVED" : "ANNOTATION_FULL_COPY";
				require(status.equals(expectedStatus), "Original disposition must match NcrnaAnnotatedCorpus.Gene.status precedence: "+id);
				final String header=id+" source="+name+" family=RF00050 status="+status;
				final Member member=new Member(genome, row);
				byLocus.put(id, member);
				byHeader.put(header, member);
				members.add(member);
				genome.observed++;
			}
		}
		require(members.size()==expected && members.size()==expectedFromRoster,
			"Original loci and roster must conserve the declared construction population");
		for(final OriginalGenome genome:originalGenomes){
			require(genome.observed==genome.expected, "Original per-genome locus census differs: "+genome.name);
		}
	}

	/** Preserve multi-interval source semantics; never silently replace them by their outer span. */
	private static int checkIntervals(final byte[] row, final int from, final int to, final int length, final int limit){
		require(from<to, "Original source intervals cannot be empty");
		long total=0;
		int previous=0;
		int count=0;
		for(int begin=from; begin<to;){
			int dash=begin, end=begin;
			while(dash<to && row[dash]>='0' && row[dash]<='9'){dash++;}
			require(dash>begin && dash-begin<=10 && dash<to && row[dash]=='-',
				"Original intervals require ordered one-based closed start-stop pairs");
			end=dash+1;
			while(end<to && row[end]>='0' && row[end]<='9'){end++;}
			require(end>dash+1 && end-dash-1<=10 && (end==to || row[end]==',' && end+1<to),
				"Original interval ends must be decimal integers separated by commas");
			final long start=Parse.parseLong(row, begin, dash), stop=Parse.parseLong(row, dash+1, end);
			require(start>previous && stop>=start && stop<=limit, "Original source intervals overlap or exceed contig bounds");
			total=Math.addExact(total, stop-start+1);
			previous=(int)stop;
			begin=end+1;
			count++;
		}
		require(total==length, "Original disjoint interval lengths must equal the extracted core length");
		return count;
	}

	private void bool(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Original provenance booleans must be literal true/false");
	}

	private void assignments(final String path){
		int count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, ASSIGNMENTS_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Original builder assignments require six fields including complete record ID");
				final Member member=byHeader.get(p.parseString(0));
				require(member!=null && member.assignment==null, "Every exact original header needs exactly one assignment");
				final String label=p.parseString(4), status=p.parseString(5);
				require(status.equals("kept") || status.equals("small") || status.equals("orphan"), "Unknown original builder status: "+status);
				require(status.equals("kept") ? !label.isEmpty() : label.isEmpty(), "Only kept original members may carry a model assignment");
				member.assignment=row;
				member.kept=status.equals("kept");
				if(member.kept){kept++;}
				count++;
			}
		}
		require(count==members.size(), "Assignments must cover all original seed inputs, including small/orphan members");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes a draft binding ledger. Every row explicitly forbids interpreting
	 * identity-only bindings as frozen training eligibility or sequence alias proof. */
	public void write(final String prefix){
		RrnaResourceIO.fresh(prefix, ".members.tsv.gz", ".bindings.tsv", ".summary.tsv");
		final ByteBuilder b=new ByteBuilder();
		try(Output rows=new Output(prefix+".members.tsv.gz"); Output bindings=new Output(prefix+".bindings.tsv");
				Output summary=new Output(prefix+".summary.tsv")){
			rows.writer.println("old_seed_input\told_consensus_member\tidentity_holdout\tsequence_alias_status\teligibility_status\t"+LOCI_HEADER
				+"\tassignment_full_record_id\tassignment_source_tid\tassignment_group\tassignment_cluster\tassignment_model\tassignment_status");
			for(final Member member:members){
				final OriginalGenome genome=member.genome;
				boolean evaluationPath=false;
				for(final Binding binding:genome.matches.bindings){evaluationPath|=binding.evaluation;}
				final String holdout=evaluationPath ? "EVALUATION_PATH" : genome.evaluationAssembly ? "EVALUATION_ASSEMBLY_CLAIM" : "NO_RECORDED_EVALUATION_IDENTITY";
				rows.writer.print(b.clear().append("true\t").append(member.kept).tab().append(holdout)
					.append("\tUNVERIFIED\tNOT_EVALUATED\t").append(member.original).tab().append(member.assignment).nl());
			}
			bindings.writer.println("original_genome\toriginal_fasta\toriginal_gff\tcohort\tinput_id\tbasis");
			long bindingCount=0;
			for(final OriginalGenome genome:originalGenomes){
				for(final Binding binding:genome.matches.bindings){
					bindings.writer.print(b.clear().append(genome.name).tab().append(genome.fasta).tab().append(genome.gff).tab()
						.append(binding.evaluation ? "evaluation" : "discovery").tab().append(binding.input)
						.append("\tLITERAL_CATALOG_PATH_NOT_SEQUENCE_PROOF\n"));
					bindingCount++;
				}
			}
			summary.writer.println("original_genomes\told_seed_members\told_consensus_members\tseed_only_members\tpath_bindings\tstatus");
			summary.writer.print(b.clear().append(originalGenomes.size()).tab().append(members.size()).tab().append(kept).tab()
				.append(members.size()-kept).tab().append(bindingCount).append("\tCONSTRUCTION_IDENTITIES_BOUND_NOT_ELIGIBILITY_OR_SEQUENCE_VALIDATION\n"));
		}
	}

	public int members(){return members.size();}
	public int kept(){return kept;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class PathBindings{
		final ArrayList<Binding> bindings=new ArrayList<Binding>();
	}
	private static final class Binding{
		Binding(final boolean evaluation_, final String input_){
			evaluation=evaluation_;
			input=input_;
		}
		final boolean evaluation;
		final String input;
		boolean evaluationAssembly;
	}
	private static final class OriginalGenome{
		OriginalGenome(final String name_, final String assembly_, final String taxid_, final String domain_, final String phylum_,
				final String gff_, final String fasta_, final int expected_, final PathBindings matches_){
			name=name_;
			assembly=assembly_;
			taxid=taxid_;
			domain=domain_;
			phylum=phylum_;
			gff=gff_;
			fasta=fasta_;
			expected=expected_;
			matches=matches_;
		}
		final String name, assembly, taxid, domain, phylum, gff, fasta;
		final int expected;
		final PathBindings matches;
		int observed;
		boolean evaluationAssembly;
	}
	private static final class Member{
		Member(final OriginalGenome genome_, final byte[] original_){
			genome=genome_;
			original=original_;
		}
		final OriginalGenome genome;
		final byte[] original;
		byte[] assignment;
		boolean kept;
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final ObjectMap<String, PathBindings> paths=new ObjectMap<String, PathBindings>(String.class, PathBindings.class);
	private final ObjectSet<String> evaluationAssemblies=new ObjectSet<String>(String.class);
	private final ObjectMap<String, OriginalGenome> byGenome=new ObjectMap<String, OriginalGenome>(String.class, OriginalGenome.class);
	private final ObjectMap<String, Member> byLocus=new ObjectMap<String, Member>(String.class, Member.class);
	private final ObjectMap<String, Member> byHeader=new ObjectMap<String, Member>(String.class, Member.class);
	private final ArrayList<OriginalGenome> originalGenomes=new ArrayList<OriginalGenome>();
	private final ArrayList<Member> members=new ArrayList<Member>();
	private long expectedFromRoster;
	private int kept;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final String GENOMES_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs";
	public static final String LOCI_HEADER="locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines";
	public static final String ASSIGNMENTS_HEADER="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus";
}
