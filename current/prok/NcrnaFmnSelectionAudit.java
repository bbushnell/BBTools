package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

import map.IntHashMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Saved original-panel membership, native standard-GA and species comparison.
 * No genome or annotation sequence is read and no new search or split is made.
 * Legacy printed-score GA_pass is deliberately ignored: actual native ! rows
 * and the completed --cut_ga command establish inclusion for retained witnesses.
 * @author Ganyu
 */
public final class NcrnaFmnSelectionAudit{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaFmnSelectionAudit.class, false);
		try{
			final Config c=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("catalog")){c.catalog=b;}
				else if(a.equals("lineages")){c.lineages=b;}
				else if(a.equals("taxonomy")){c.taxonomy=b;}
				else if(a.equals("hits")){c.hits=b;}
				else if(a.equals("expectedgenomes")){c.pool=Parse.parseIntKMG(b);}
				else if(a.equals("expectedevaluation")){c.evaluation=Parse.parseIntKMG(b);}
				else if(a.equals("expectedoverlap")){c.overlap=Parse.parseIntKMG(b);}
				else if(a.equals("expectedlineages")){c.labels=Parse.parseIntKMG(b);}
				else if(a.equals("expectedhits")){c.hitRows=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown selection audit option: "+arg);}
			}
			c.roster=parser.in1;
			c.output=parser.out1;
			require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Audit complete saved populations without sampling or overwrite");
			run(c);
			pp.outstream.println("NCRNA_FMN_SELECTION_AUDIT_COMPLETE_NOT_SPLIT");
		}finally{Shared.closeStream(pp.outstream);}
	}
	public static void run(final Config c){
		require(c!=null && c.roster!=null && c.catalog!=null && c.lineages!=null && c.taxonomy!=null
			&& c.hits!=null && c.output!=null && c.pool>0 && c.evaluation>0 && c.overlap>=0 && c.labels>0,
			"Declare original pool, earlier scan catalog, saved lineages, taxonomy, original standard hits and expected populations");
		RrnaResourceIO.fresh(c.output, ".inputs.tsv", ".species.tsv", ".summary.tsv", ".witnesses.tsv");
		final NcrnaFmnSelectionAudit x=new NcrnaFmnSelectionAudit();
		x.roster(c.roster, c.pool);
		x.catalog(c.catalog, c.evaluation, c.overlap);
		x.lineages(c.lineages);
		require(x.items.size()==c.labels, "The saved lineage population must be complete");
		x.taxonomy(c.taxonomy);
		x.hits(c.hits, c.output, c.hitRows);
		x.write(c);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void roster(final String path, final int expected){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GENOME_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				require(pool.add(text(7)), "Original PGAP paths must be unique");
			}
		}
		require(pool.size()==expected, "The original PGAP roster is incomplete");
	}
	private void catalog(final String path, final int expected, final int overlap){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, CENSUS_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 5);
				final int id=p.parseInt(0);
				final String fasta=text(1);
				require(id>0 && !catalogIds.contains(id) && catalogPaths.add(fasta), "Earlier input IDs and paths must be unique");
				catalogIds.put(id, earlier.size());
				earlier.add(fasta);
				if(pool.contains(fasta)){earlierPool++;}
			}
		}
		require(earlier.size()==expected && earlierPool==overlap, "The earlier PGAP intersection must reproduce the saved path subtraction");
	}
	private void lineages(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, LINEAGE_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 11);
				final Item x=new Item();
				x.cohort=text(0);
				x.id=p.parseInt(1);
				x.path=text(2);
				x.group=p.parseInt(3);
				x.assembly=text(4);
				x.genus=text(8);
				x.cores=p.parseInt(10);
				x.inPool=pool.contains(x.path);
				require(x.id>0 && x.group>=0 && x.cores>=0 && p.termEquals("GTDB_ACCESSION_MATCH", 9)
					&& (x.cohort.equals("evaluation") || x.cohort.equals("discovery")), "Saved lineage identity must be complete");
				if(x.cohort.equals("evaluation")){
					final int index=catalogIds.get(x.id);
					require(index>=0 && earlier.get(index).equals(x.path) && !evaluationItems.contains(x.id), "Evaluation labels must bind original scan ID and path");
					evaluationItems.put(x.id, items.size());
				}else{require(x.inPool && !catalogPaths.contains(x.path), "Discovery rows must belong to the literal PGAP remainder");}
				require(itemKeys.add(x.cohort+'\t'+x.id), "Lineage rows must be unique");
				items.add(x);
			}
		}
	}
	private void taxonomy(final String path){
		final ObjectMap<NcrnaReferenceGeometry.Target, String> species=new ObjectMap<NcrnaReferenceGeometry.Target, String>(NcrnaReferenceGeometry.Target.class, String.class);
		final NcrnaReferenceGeometry.Target key=new NcrnaReferenceGeometry.Target();
		for(final Item x:items){final byte[] a=utf8(x.assembly); key.set(a, 0, a.length); species.put(key.copy(), null);}
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 2);
				p.setBounds(0);
				int start=p.a();
				if(p.b()-start>3 && row[start+2]=='_' && ((row[start]=='R' && row[start+1]=='S') || (row[start]=='G' && row[start+1]=='B'))){start+=3;}
				key.set(row, start, p.b());
				if(!species.contains(key)){continue;}
				require(species.get(key)==null, "Duplicate GTDB accession would make the species join ambiguous");
				final String[] ranks=text(1).split(";", -1);
				require(ranks.length==7 && ranks[1].equals("p__Bacteroidota") && ranks[5].startsWith("g__") && ranks[6].startsWith("s__"), "Saved Bacteroidota labels must match original GTDB taxonomy");
				species.put(key.copy(), ranks[5]+'\t'+ranks[6]);
			}
		}
		for(final Item x:items){
			final byte[] a=utf8(x.assembly);
			key.set(a, 0, a.length);
			final String value=species.get(key);
			require(value!=null && value.startsWith("g__"+x.genus+'\t'), "Original GTDB genus must reproduce every saved lineage label");
			x.species=value.substring(value.indexOf('\t')+1);
			if(x.species.equals("s__")){x.species="UNRESOLVED_SPECIES:"+x.assembly;}
		}
	}
	private void hits(final String path, final String output, final int expected){
		final ObjectMap<String, Table> byTable=new ObjectMap<String, Table>(String.class, Table.class);
		final ArrayList<String> tables=new ArrayList<String>();
		long total=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.HIT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 21);
				total++;
				require(catalogIds.contains(p.parseInt(0)), "Standard-scan hit must bind an earlier input");
				if(!p.termEquals("RF00050", 7)){continue;}
				final int index=evaluationItems.get(p.parseInt(0));
				if(index<0){continue;}
				final String table=text(1);
				require(table.startsWith("/"), "Native witnesses require literal saved paths");
				Table entry=byTable.get(table);
				if(entry==null){
					entry=new Table();
					byTable.put(table, entry);
					tables.add(table);
				}
				entry.rows.add(new Witness(items.get(index), p.parseLong(2), text(3), p.parseInt(4), p.parseInt(5), text(6)));
				// Field17 is historical printed-score GA_pass; never use it as native inclusion.
			}
		}
		require(total==expected, "Original standard hit rows must conserve the accepted complete evidence count");
		Collections.sort(tables);
		final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();
		try(Output out=new Output(output+".witnesses.tsv")){
			out.writer.println("input_id\tnative_table\tnative_line\ttarget\tstart1\tstop1\tstrand\tnative_included");
			for(final String table:tables){
				final ArrayList<Witness> list=byTable.get(table).rows;
				Collections.sort(list, Comparator.comparingLong(w->w.line));
				int next=0;
				long number=0;
				boolean ga=false, done=false;
				try(Input in=new Input(table)){
					for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
						number++;
						if(row.length>0 && row[0]=='#'){
							final String comment=new String(row, StandardCharsets.US_ASCII);
							ga|=comment.startsWith("# Option settings:") && comment.contains(" --cut_ga ");
							done|=comment.equals("# [ok]");
						}
						if(next==list.size() || number!=list.get(next).line){continue;}
						final Witness w=list.get(next++);
						f.set(row);
						require(f.eq(2, "FMN") && f.eq(3, "RF00050") && f.eq(4, "cm") && f.eq(16, "!")
							&& f.eq(0, w.target) && f.eq(9, w.strand) && Math.min(f.number(7), f.number(8))==w.start
							&& Math.max(f.number(7), f.number(8))==w.stop, "Original standard-GA physical row must match input target/strand/bounds and native inclusion");
						w.item.ga++;
						out.writer.print(b.clear().append(w.item.id).tab().append(utf8(table)).tab().append(number).tab().append(utf8(w.target))
							.tab().append(w.start).tab().append(w.stop).tab().append(w.strand).append("\ttrue\n"));
					}
				}
				require(next==list.size() && ga && done, "Every earlier native witness requires complete standard-GA output; duplicate or missing physical lines fail");
			}
		}
	}
	private void write(final Config c){
		final ObjectMap<String, Species> species=new ObjectMap<String, Species>(String.class, Species.class);
		final ArrayList<Species> ordered=new ArrayList<Species>();
		final ObjectMap<String, String> groupSpecies=new ObjectMap<String, String>(String.class, String.class);
		final long[] counts=new long[4];
		try(Output inputs=new Output(c.output+".inputs.tsv"); Output summary=new Output(c.output+".summary.tsv"); Output speciesOut=new Output(c.output+".species.tsv")){
			inputs.writer.println("cohort\tinput_id\tfasta\tassembly_group\tassembly\tgenus\tspecies\tin_original_PGAP\tCM_admitted_cores\tearlier_standard_GA_FMN_events");
			for(final Item x:items){
				final boolean eval=x.cohort.equals("evaluation");
				inputs.writer.print(b.clear().append(x.cohort).tab().append(x.id).tab().append(utf8(x.path)).tab().append(x.group)
					.tab().append(x.assembly).tab().append(utf8(x.genus)).tab().append(utf8(x.species)).tab().append(x.inPool).tab().append(x.cores).tab());
				if(eval){inputs.writer.print(b.clear().append(x.ga).nl());}else{inputs.writer.println("NA");}
				if(eval && x.inPool){counts[0]++; if(x.ga>0){counts[1]++;}}
				if(!x.genus.equals("Chryseobacterium")){continue;}
				final String groupKey=x.cohort+'\t'+x.group;
				final String previous=groupSpecies.put(groupKey, x.species);
				require(previous==null || previous.equals(x.species), "One resolved Chryseobacterium group has conflicting GTDB species");
				if(eval && x.inPool){counts[2]++; if(x.ga>0){counts[3]++;}}
				Species s=species.get(x.species);
				if(s==null){s=new Species(x.species); species.put(x.species, s); ordered.add(s);}
				final IntHashMap set=eval ? x.cores>0 ? s.positive : s.negative : s.discovery;
				set.put(x.group, 1);
			}
			Collections.sort(ordered, Comparator.comparing(s->s.name));
			speciesOut.writer.println("GTDB_species\tevaluation_clean_core_positive_groups\tevaluation_without_clean_core_groups\tdiscovery_groups\tnegative_evaluation_and_discovery");
			int shared=0, negativeShared=0;
			for(final Species s:ordered){
				for(final int group:s.positive.keys()){if(group>=0){s.negative.remove(group);}}
				final boolean overlap=!s.name.startsWith("UNRESOLVED_SPECIES:") && s.negative.size()>0 && s.discovery.size()>0;
				if(overlap){shared++; negativeShared+=s.negative.size();}
				speciesOut.writer.print(b.clear().append(utf8(s.name)).tab().append(s.positive.size()).tab().append(s.negative.size())
					.tab().append(s.discovery.size()).tab().append(overlap).nl());
			}
			summary.writer.println("original_PGAP_inputs\tearlier_scan_inputs\tearlier_PGAP_inputs\tremaining_PGAP_inputs\tearlier_PGAP_Bacteroidota\tBacteroidota_standard_GA_positive_inputs\tearlier_PGAP_Chryseobacterium\tChryseobacterium_standard_GA_positive_inputs\tnegative_evaluation_species_shared_with_discovery\tnegative_evaluation_groups_in_shared_species\tstatus");
			b.clear().append(pool.size()).tab().append(earlier.size()).tab().append(earlierPool).tab().append(pool.size()-earlierPool);
			for(final long value:counts){b.tab().append(value);}
			summary.writer.print(b.tab().append(shared).tab().append(negativeShared).append("\tSAVED_SELECTION_DIAGNOSTIC_NOT_SPLIT\n"));
		}
	}
	private void parse(final byte[] row, final int fields){p.set(row); require(p.terms()==fields, "Saved selection metadata schema changed");}
	private String text(final int field){return p.parseString(field);}
	private static byte[] utf8(final String s){return s.getBytes(StandardCharsets.UTF_8);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String roster, catalog, lineages, taxonomy, hits, output;
		public int pool=39979, evaluation=11769, overlap=10666, labels=4395, hitRows=932681;
	}
	private static final class Item{
		String cohort, path, assembly, genus, species;
		int id, group, cores, ga;
		boolean inPool;
	}
	private static final class Witness{
		Witness(final Item x, final long n, final String t, final int a, final int z, final String s){
			require(n>0 && a>0 && z>=a && (s.equals("+") || s.equals("-")), "A native witness needs physical line and valid original coordinates");
			item=x;
			line=n;
			target=t;
			start=a;
			stop=z;
			strand=s;
		}
		final Item item;
		final long line;
		final String target, strand;
		final int start, stop;
	}
	private static final class Species{
		Species(final String n){name=n;}
		final String name;
		final IntHashMap positive=new IntHashMap(), negative=new IntHashMap(), discovery=new IntHashMap();
	}
	private static final class Table{
		final ArrayList<Witness> rows=new ArrayList<Witness>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final ObjectSet<String> pool=new ObjectSet<String>(String.class), catalogPaths=new ObjectSet<String>(String.class), itemKeys=new ObjectSet<String>(String.class);
	private final IntHashMap catalogIds=new IntHashMap(), evaluationItems=new IntHashMap();
	private final ArrayList<String> earlier=new ArrayList<String>();
	private final ArrayList<Item> items=new ArrayList<Item>();
	private int earlierPool;
	private static final String CENSUS_HEADER="input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records";
	private static final String LINEAGE_HEADER="cohort\tinput_id\tfasta\tassembly_group\tassembly\tclass\torder\tfamily\tgenus\ttaxonomy_status\tCM_admitted_cores";
}
