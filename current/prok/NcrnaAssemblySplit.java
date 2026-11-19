package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import java.util.Random;

import map.IntHashMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** One hit-blind genus-stratified assembly split shared by all RNA families.
 * Original roster locus counts and historical evaluation membership are not read
 * into selection. Rare strata train wholly; ordinary strata hold out floor(n/5).
 * @author Ganyu
 */
public final class NcrnaAssemblySplit{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaAssemblySplit.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String groups=null, taxonomy=null, coassign=null;
			long seed=-1;
			int expected=39979;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("groups")){groups=b;}
				else if(a.equals("taxonomy")){taxonomy=b;}
				else if(a.equals("coassign")){coassign=b;}
				else if(a.equals("seed")){seed=Long.parseLong(b);}
				else if(a.equals("expectedgenomes")){expected=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown assembly split option: "+arg);}
			}
			require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Split complete populations without sampling, overwrite or paired inputs");
			run(parser.in1, groups, taxonomy, parser.out1, seed, expected, coassign);
			pp.outstream.println("NCRNA_ASSEMBLY_SPLIT_COMPLETE_HIT_BLIND_NOT_MEMBER_ADMISSION");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void run(final String roster, final String groups, final String taxonomy,
			final String output, final long seed, final int expected){
		run(roster, groups, taxonomy, output, seed, expected, null);
	}
	public static void run(final String roster, final String groups, final String taxonomy,
			final String output, final long seed, final int expected, final String coassign){
		require(roster!=null && groups!=null && taxonomy!=null && output!=null && seed>=0 && expected>0,
			"Declare the complete original roster, resolved groups, taxonomy files, fresh output and explicit seed");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Assembly splitting requires assertions for group conservation");
		RrnaResourceIO.fresh(output, ".inputs.tsv.gz", ".groups.tsv", ".strata.tsv", ".summary.tsv", ".coassigned.tsv");
		final NcrnaAssemblySplit x=new NcrnaAssemblySplit();
		x.roster(roster, expected);
		if(coassign!=null){x.coassignment(coassign);}
		x.groups(groups);
		for(final String file:taxonomy.split(",")){x.taxonomy(file);}
		x.assign(seed);
		x.write(output, seed);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void roster(final String path, final int expected){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GENOME_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				final Genome g=new Genome(text(2), text(7));
				require(!g.accession.isEmpty() && g.path.startsWith("/") && byPath.get(g.path)==null,
					"Every original genome has one nonempty accession and unique literal path");
				byPath.put(g.path, g);
				genomes.add(g);
				final byte[] accession=utf8(g.accession);
				probe.set(accession, 0, accession.length);
				require(!byAccession.contains(probe), "Original assembly accessions must be unique");
				byAccession.put(probe.copy(), g);
				// Deliberately do not parse source-locus count, taxonomy guesses, or annotations.
			}
		}
		require(genomes.size()==expected, "Original pool size differs from its frozen manifest");
	}

	private void groups(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GROUP_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 14);
				final int id=p.parseInt(10);
				final Genome g=byPath.get(text(2));
				final boolean exception=id==coassignedGroup || id==partnerGroup;
				if(exception){
					pairRows.add(new PairRow(id, text(2), g!=null));
					if(id==partnerGroup){partnerSeen=true;}
				}
				if(g==null){continue;}// Other historical paths are outside the declared PGAP pool.
				final boolean clean=p.termEquals("false", 12) && !text(13).contains("UNRESOLVED") && !text(13).contains("PENDING");
				final boolean allowed=exception && p.termEquals("true", 12)
					&& (p.termEquals("ASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED", 13)
					|| p.termEquals("EXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS", 13) || p.termEquals("EVALUATION", 13));
				require(g.group==null && id>=0 && (clean || allowed),
					"Resolve every pool assembly before splitting; duplicate/unresolved source: "+g.path);
				coassignedUnresolved|=allowed;
				int index=groupIndex.get(id);
				if(index<0){
					index=allGroups.size();
					groupIndex.put(id, index);
					allGroups.add(new Group(id));
				}
				g.group=allGroups.get(index);
				g.group.members.add(g);
			}
		}
		for(final Genome g:genomes){require(g.group!=null, "Pool genome lacks a resolved assembly group: "+g.path);}
		if(coassignedGroup>=0){
			require(groupIndex.get(coassignedGroup)>=0 && coassignedUnresolved && partnerSeen, "The explicit unresolved co-assignment must bind both saved groups");
			coassigned=allGroups.get(groupIndex.get(coassignedGroup));
			final int partnerIndex=groupIndex.get(partnerGroup);
			partner=partnerIndex<0 ? null : allGroups.get(partnerIndex);
		}
	}

	/** Narrow reviewed exception: preserve both groups and record their same side.
	 * Outside-pool aliases inherit that side without becoming additional quota units. */
	private void coassignment(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "pool_group\tpartner_group\treason", path);
			final byte[] row=in.reader.nextLine();
			require(row!=null, "An explicit co-assignment file needs its reviewed pair");
			parse(row, 3);
			coassignedGroup=p.parseInt(0);
			partnerGroup=p.parseInt(1);
			require(coassignedGroup>=0 && partnerGroup>=0 && coassignedGroup!=partnerGroup
				&& p.termEquals("unresolved, co-assigned to prevent leakage", 2) && in.reader.nextLine()==null,
				"Only one expressly reviewed pair is supported; no identity or further unresolved exceptions are inferred");
		}
	}

	private void taxonomy(final String path){
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 2);
				p.setBounds(0);
				int start=p.a();
				if(p.b()-start>3 && row[start+2]=='_' && ((row[start]=='R' && row[start+1]=='S') || (row[start]=='G' && row[start+1]=='B'))){start+=3;}
				probe.set(row, start, p.b());
				final Genome g=byAccession.get(probe);
				if(g==null){continue;}
				require(g.stratum==null, "One original accession cannot receive two taxonomy assignments");
				g.stratum=stratum(text(1));
			}
		}
	}

	/** Full lineage prefixes prevent equal genus/family text in different lineages colliding. */
	public static String stratum(final String taxonomy){
		final String[] ranks=taxonomy.split(";", -1);
		final String[] prefixes={"d__", "p__", "c__", "o__", "f__", "g__", "s__"};
		require(ranks.length==7, "GTDB taxonomy requires exactly seven ranked fields");
		for(int i=0; i<ranks.length; i++){require(ranks[i].startsWith(prefixes[i]), "GTDB rank order changed");}
		final int last=ranks[5].length()>3 ? 5 : 4;
		require(ranks[last].length()>3, "Missing genus and family require explicit taxonomy adjudication");
		final StringBuilder out=new StringBuilder(last==5 ? "GENUS:" : "FAMILY:");
		for(int i=0; i<=last; i++){
			if(i>0){out.append(';');}
			out.append(ranks[i]);
		}
		return out.toString();
	}

	private void assign(final long seed){
		final ObjectMap<String, Stratum> byStratum=new ObjectMap<String, Stratum>(String.class, Stratum.class);
		for(final Group group:allGroups){
			for(final Genome g:group.members){
				require(g.stratum!=null, "Original accession has no GTDB genus/family: "+g.accession);
				if(group.stratum==null){group.stratum=g.stratum;}
				require(group.stratum.equals(g.stratum), "One assembly group has conflicting GTDB strata: "+group.id);
			}
			Stratum s=byStratum.get(group.stratum);
			if(s==null){
				s=new Stratum(group.stratum);
				byStratum.put(s.name, s);
				strata.add(s);
			}
			s.groups.add(group);
		}
		Collections.sort(strata, Comparator.comparing(s->s.name));
		final Random random=new Random(seed);
		for(final Stratum s:strata){
			Collections.sort(s.groups, Comparator.comparingInt(g->g.id));
			final boolean rare=s.groups.size()<5;
			final int held=rare ? 0 : s.groups.size()/5;
			Collections.shuffle(s.groups, random);
			for(int i=0; i<s.groups.size(); i++){
				final Group g=s.groups.get(i);
				g.assignment=rare ? "COVERAGE_CRITICAL_TRAIN" : i<held ? "HELD_OUT" : "TRAIN";
				if(i<held){s.held++;}
				s.inputs+=g.members.size();
			}
			require(s.held==held, "Each stratum must conserve its exact declared heldout group count");
		}
		if(partner!=null){
			require(partner.stratum.equals(coassigned.stratum) && partner.assignment.equals(coassigned.assignment),
				"The reviewed pair must share a stratum and side without altering group quotas; a changed population needs an explicit new disposition");
		}
		Collections.sort(allGroups, Comparator.comparingInt(g->g.id));
		Collections.sort(genomes, Comparator.comparing(g->g.accession));
	}

	private void write(final String prefix, final long seed){
		final ByteBuilder b=new ByteBuilder();
		long held=0, critical=0, train=0;
		try(Output input=new Output(prefix+".inputs.tsv.gz"); Output groupOut=new Output(prefix+".groups.tsv");
				Output strataOut=new Output(prefix+".strata.tsv"); Output summary=new Output(prefix+".summary.tsv");
				Output paired=new Output(prefix+".coassigned.tsv")){
			paired.writer.println("pool_group\tpartner_group\tinput_group\tfasta\tin_PGAP_pool\tassignment\treason");
			Collections.sort(pairRows, Comparator.comparing(row->row.path));
			for(final PairRow row:pairRows){
				paired.writer.print(b.clear().append(coassignedGroup).tab().append(partnerGroup).tab().append(row.group)
					.tab().append(utf8(row.path)).tab().append(row.inPool)
					.tab().append(coassigned.assignment).append("\tunresolved, co-assigned to prevent leakage\n"));
			}
			input.writer.println("assembly\tfasta\tassembly_group\tstratum\tassignment");
			for(final Genome g:genomes){
				input.writer.print(b.clear().append(utf8(g.accession)).tab().append(utf8(g.path)).tab().append(g.group.id)
					.tab().append(utf8(g.stratum)).tab().append(g.group.assignment).nl());
			}
			groupOut.writer.println("assembly_group\tstratum\tpool_inputs\tassignment");
			for(final Group g:allGroups){
				require(g.assignment!=null, "Every resolved pool group must have exactly one assignment");
				groupOut.writer.print(b.clear().append(g.id).tab().append(utf8(g.stratum)).tab().append(g.members.size()).tab().append(g.assignment).nl());
				if(g.assignment.equals("HELD_OUT")){held++;}
				else if(g.assignment.equals("COVERAGE_CRITICAL_TRAIN")){critical++;}
				else{train++;}
			}
			strataOut.writer.println("stratum\tpool_inputs\tassembly_groups\theldout_groups\ttraining_groups\tcoverage_critical");
			for(final Stratum s:strata){
				strataOut.writer.print(b.clear().append(utf8(s.name)).tab().append(s.inputs).tab().append(s.groups.size())
					.tab().append(s.held).tab().append(s.groups.size()-s.held).tab().append(s.groups.size()<5).nl());
			}
			require(held+critical+train==allGroups.size(), "Mutually exclusive assignments must conserve every pool group");
			summary.writer.println("pool_inputs\tassembly_groups\tstrata\theldout_groups\tordinary_train_groups\tcoverage_critical_groups\tseed\trule\tstatus");
			summary.writer.print(b.clear().append(genomes.size()).tab().append(allGroups.size()).tab().append(strata.size())
				.tab().append(held).tab().append(train).tab().append(critical).tab().append(seed)
				.append("\tSORTED_STRATA_SORTED_GROUP_IDS_JAVA_RANDOM_SHUFFLE_FLOOR_N_OVER5\tHIT_BLIND_GROUP_SPLIT_NOT_MEMBER_ADMISSION\n"));
		}
	}

	private void parse(final byte[] row, final int fields){p.set(row); require(p.terms()==fields, "Saved metadata schema changed");}
	private String text(final int field){return p.parseString(field);}
	private static byte[] utf8(final String s){return s.getBytes(StandardCharsets.UTF_8);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome{
		Genome(final String a, final String p){accession=a; path=p;}
		final String accession, path;
		String stratum;
		Group group;
	}
	private static final class Group{
		Group(final int value){id=value;}
		final int id;
		final ArrayList<Genome> members=new ArrayList<Genome>();
		String stratum, assignment;
	}
	private static final class Stratum{
		Stratum(final String value){name=value;}
		final String name;
		final ArrayList<Group> groups=new ArrayList<Group>();
		int held, inputs;
	}
	private static final class PairRow{
		PairRow(final int id, final String fasta, final boolean pool){
			group=id;
			path=fasta;
			inPool=pool;
		}
		final int group;
		final String path;
		final boolean inPool;
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final NcrnaReferenceGeometry.Target probe=new NcrnaReferenceGeometry.Target();
	private final ObjectMap<String, Genome> byPath=new ObjectMap<String, Genome>(String.class, Genome.class);
	private final ObjectMap<NcrnaReferenceGeometry.Target, Genome> byAccession=
		new ObjectMap<NcrnaReferenceGeometry.Target, Genome>(NcrnaReferenceGeometry.Target.class, Genome.class);
	private final IntHashMap groupIndex=new IntHashMap();
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Group> allGroups=new ArrayList<Group>();
	private final ArrayList<Stratum> strata=new ArrayList<Stratum>();
	private int coassignedGroup=-1, partnerGroup=-1;
	private Group coassigned, partner;
	private boolean coassignedUnresolved, partnerSeen;
	private final ArrayList<PairRow> pairRows=new ArrayList<PairRow>();
}
