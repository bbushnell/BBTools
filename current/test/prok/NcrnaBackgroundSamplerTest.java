package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CobalaminShredSampler;
import prok.Euk18sShredSet;
import prok.NcrnaCmLocusAudit;
import prok.NcrnaCorpusSample;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.RrnaResourceIO;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Literal GC/interval oracles and native panel/sampling CLI controls.
 * All sequence and table I/O uses the same BBTools owners as the public tools.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBackgroundSamplerTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaBackgroundSamplerTest.class, false);
		try{
			final Parser parser=new Parser();
			for(String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			check(parser.out1!=null, "out= must name a fresh fixture directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be new: "+root);
			gcCases();
			intervalCases();
			final Fixture fixture=new Fixture(root);
			fixture.write(false);
			panel(fixture, "panel", 5, fixture.genomes, fixture.groups);
			checkPanel(fixture, "panel");
			sample(fixture, "sample", fixture.ids, 10, 10, fixture.lengths, "stratify=phylumgc", "gcbins=10", "nonoverlap=t");
			inspect(fixture, "sample", 10, true);
			sample(fixture, "selected_panel", path(root, "panel.ids.txt"), 5, 20, fixture.lengths,
				"stratify=phylumgc", "gcbins=3", "nonoverlap=t");
			inspect(fixture, "selected_panel", 20, true);
			reverse(fixture.genomes, path(root, "genomes.reverse.tsv"), true);
			reverse(fixture.groups, path(root, "groups.reverse.tsv"), true);
			reverse(fixture.ids, path(root, "ids.reverse.txt"), false);
			fixture.writeSequences(true);
			panel(fixture, "panel_reverse", 5, path(root, "genomes.reverse.tsv"), path(root, "groups.reverse.tsv"));
			for(String suffix:new String[]{".ids.txt", ".genomes.tsv", ".groups.tsv", ".summary.tsv"}){
				equalLines(path(root, "panel"+suffix), path(root, "panel_reverse"+suffix));
			}
			sample(fixture, "reverse", path(root, "ids.reverse.txt"), 10, 10, fixture.lengths,
				"stratify=phylumgc", "gcbins=10", "nonoverlap=t", "genomes="+path(root, "genomes.reverse.tsv"));
			for(String suffix:new String[]{".cores.fa.gz", ".contexts.fa.gz", ".shreds.tsv", ".genomes.tsv", ".lengths.tsv",
				".summary.tsv", ".schedule.tsv", ".window_metadata.tsv", ".strata.tsv"}){
				equalLines(path(root, "sample"+suffix), path(root, "reverse"+suffix));
			}
			reject(()->panel(fixture, "too_small", 2, fixture.genomes, fixture.groups), "cannot cover every promised phylum");
			check(!new File(root, "too_small.ids.txt").exists(), "Impossible panel quotas fail before output");
			writeLengths(path(root, "long.fa"), new int[]{600, 600, 600, 600});
			reject(()->sample(fixture, "exhausted", fixture.ids, 10, 10, path(root, "long.fa"),
				"stratify=phylumgc", "nonoverlap=t"), "Fixed schedule exhausted:");
			check(countRows(path(root, "exhausted.schedule.tsv"))==10, "Exhaustion preserves every requested source/stratum/length");
			final int[] partial={0};
			RrnaResourceIO.stream(path(root, "exhausted.cores.fa.gz"), read->partial[0]++);
			check(partial[0]>0 && partial[0]<10, "Partial FASTA closes cleanly and never pretends all exhausted targets were emitted");
			reject(()->sample(fixture, "bad_bins", fixture.ids, 10, 10, fixture.lengths, "stratify=phylumgc", "gcbins=0"), "gcbins requires");
			check(!new File(root, "bad_bins.schedule.tsv").exists(), "Invalid modes fail before output");
			write(path(root, "duplicate.ids"), "g0\ng0\n");
			reject(()->sample(fixture, "duplicate", path(root, "duplicate.ids"), 2, 10, fixture.lengths), "Panel IDs must be unique");
			reject(()->sample(fixture, "mixed_mode", fixture.ids, 10, 10, fixture.lengths, "groups="+fixture.groups),
				"groups/panelcount belong to mode=panel");
			copyWithoutLast(fixture.genomes, path(root, "partial_group.tsv"));
			reject(()->CobalaminShredSampler.main(new String[]{"mode=panel", "genomes="+path(root, "partial_group.tsv"),
				"groups="+fixture.groups, "outprefix="+path(root, "partial_group"), "expectedgenomes=10", "seed=47"}),
				"Source groups must enter or leave the candidate roster whole");
			final byte[] undefined=new byte[1000];
			Arrays.fill(undefined, (byte)'N');
			writeGenome(path(root, "g9.fa"), undefined, false);
			reject(()->sample(fixture, "undefined", fixture.ids, 10, 10, fixture.lengths, "stratify=phylumgc"), "without defined bases: g9");
			pp.outstream.println("NCRNA_BACKGROUND_SAMPLER_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Exact integer boundary cases would detect rounding before bin assignment. */
	private static void gcCases(){
		check(CobalaminShredSampler.gcBin(0, 100, 10)==0, "Zero GC belongs in the first bin");
		check(CobalaminShredSampler.gcBin(39, 100, 10)==3, "Below 40% remains below the boundary");
		check(CobalaminShredSampler.gcBin(40, 100, 10)==4, "A half-open 40% boundary enters bin four");
		check(CobalaminShredSampler.gcBin(100, 100, 10)==9, "100% belongs in the last bin");
		check(CobalaminShredSampler.gcBin(33, 99, 3)==1, "Configurable thirds use exact count arithmetic");
		check(CobalaminShredSampler.gcBin(4000000000L, 10000000000L, 10)==4, "Whole-genome counters exceed signed-int range");
		reject(()->CobalaminShredSampler.gcBin(0, 0, 10), "GC bins require");
	}

	/** A Boolean oracle checks every inclusive interval after adjacency and bridge updates. */
	private static void intervalCases(){
		final Euk18sShredSet.Contig mask=new Euk18sShredSet.Contig("literal", 40, 0, 1, 1);
		mask.mask.add(new int[]{5, 8});
		mask.mask.add(new int[]{20, 22});
		mask.prepareMask(0);
		final boolean[] excluded=new boolean[41];
		Arrays.fill(excluded, 5, 9, true);
		Arrays.fill(excluded, 20, 23, true);
		for(int[] interval:new int[][]{{9, 12}, {15, 19}, {13, 14}, {1, 2}, {40, 40}, {3, 4}, {23, 39}}){
			CobalaminShredSampler.addInterval(mask, interval[0], interval[1]);
			Arrays.fill(excluded, interval[0], interval[1]+1, true);
			for(int lo=1; lo<=40; lo++){
				boolean expected=false;
				for(int hi=lo; hi<=40; hi++){
					expected|=excluded[hi];
					check(mask.overlaps(lo, hi)==expected, "Primitive interval union must match the independent Boolean oracle");
				}
			}
		}
		reject(()->CobalaminShredSampler.addInterval(mask, 0, 1), "inside its contig");
	}

	private static void panel(final Fixture fixture, final String name, final int cap, final String genomes, final String groups){
		CobalaminShredSampler.main(new String[]{"mode=panel", "genomes="+genomes, "groups="+groups,
			"outprefix="+path(fixture.root, name), "expectedgenomes=11", "panelcount="+cap, "seed=47"});
	}

	/** The fixture has eight P groups, one Q and one R; a five-group panel must retain 3/1/1. */
	private static void checkPanel(final Fixture fixture, final String prefix){
		final ObjectIntMap<String> phyla=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		int count=0;
		try(Input input=new Input(path(fixture.root, prefix+".genomes.tsv"))){
			RrnaResourceIO.header(input.reader, NcrnaCorpusSample.GENOME_HEADER, "panel");
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				check(!p.termEquals("g10", 1), "The lexical representative of the g0/g10 source group is g0");
				final String phylum=p.parseString(5);
				phyla.put(phylum, phyla.contains(phylum) ? phyla.get(phylum)+1 : 1);
				count++;
			}
		}
		check(count==5 && phyla.get("P")==3 && phyla.get("Q")==1 && phyla.get("R")==1,
			"Reserved singleton phyla and residual-capacity weighting must realize literal 3/1/1 quotas");
	}

	private static void sample(final Fixture fixture, final String name, final String ids, final int sources, final int count,
			final String lengths, final String... extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("genomes="+fixture.genomes, "panelids="+ids,
			"lengths="+lengths, "loci="+fixture.loci, "outprefix="+path(fixture.root, name), "arm=panel", "count="+count,
			"genomecount="+sources, "expectedgenomes=11", "expectedpanel="+sources, "expectedlengths=4", "expectedloci=11",
			"families=RF00059", "seed=101", "flank=20"));
		args.addAll(Arrays.asList(extra));
		CobalaminShredSampler.main(args.toArray(new String[args.size()]));
	}

	/** Checks real extracted bytes, target lengths, masks, nonoverlap and quota accounting. */
	private static void inspect(final Fixture fixture, final String prefix, final int expected, final boolean disjoint){
		final ObjectIntMap<String> schedule=new ObjectIntMap<String>(String.class);
		final ObjectMap<String, String> sources=new ObjectMap<String, String>(String.class, String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path(fixture.root, prefix+".schedule.tsv"))){
			input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				check(!schedule.contains(p.parseString(0)), "Fixed schedule IDs must be unique");
				schedule.put(p.parseString(0), p.parseInt(3));
				sources.put(p.parseString(0), p.parseString(1));
			}
		}
		final ObjectMap<String, byte[]> cores=new ObjectMap<String, byte[]>(String.class, byte[].class);
		RrnaResourceIO.stream(path(fixture.root, prefix+".cores.fa.gz"), read->{
			check(cores.put(read.id, read.bases)==null, "Core IDs must be unique");
		});
		final ObjectMap<String, byte[]> contexts=new ObjectMap<String, byte[]>(String.class, byte[].class);
		RrnaResourceIO.stream(path(fixture.root, prefix+".contexts.fa.gz"), read->{
			check(contexts.put(read.id, read.bases)==null, "Context IDs must be unique");
		});
		final ObjectMap<String, IntList> intervals=new ObjectMap<String, IntList>(String.class, IntList.class);
		int seen=0;
		try(Input input=new Input(path(fixture.root, prefix+".shreds.tsv"))){
			RrnaResourceIO.header(input.reader, CobalaminShredSampler.HEADER, "shreds");
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				final String id=p.parseString(0), genome=p.parseString(2), contig=p.parseString(5), key=genome+"/"+contig;
				final int lo=p.parseInt(6), hi=p.parseInt(7), length=p.parseInt(12);
				check(schedule.contains(id) && schedule.get(id)==length && length==hi-lo+1, "No sampled length may be redrawn or shortened");
				check(genome.equals(sources.get(id)), "Every scheduled window must retain its assigned source");
				check(lo>=1 && hi<=1000 && (hi<201 || lo>205), "Cores stay in the real contig and exclude ambiguous bases");
				check(!contig.equals("a") || hi<50 || lo>70, "Known TPP coordinates remain excluded");
				check(Arrays.equals(cores.get(id), Arrays.copyOfRange(fixture.bases.get(genome), lo-1, hi)), "Output bases must match literal source extraction");
				final int left=Math.max(1, lo-20), right=Math.min(1000, hi+20);
				check(p.parseInt(8)==left && p.parseInt(9)==right && p.parseInt(10)==lo-left+1 && p.parseInt(11)==hi-left+1,
					"Contexts must clip at the contig boundary while preserving core coordinates");
				check(Arrays.equals(contexts.get(id), Arrays.copyOfRange(fixture.bases.get(genome), left-1, right)),
					"Context bases must match their recorded source coordinates");
				IntList prior=intervals.get(key);
				if(prior==null){
					prior=new IntList();
					intervals.put(key, prior);
				}
				if(disjoint){
					for(int i=0; i<prior.size; i+=2){check(hi<prior.get(i) || lo>prior.get(i+1), "Selected cores may not overlap across different scheduled lengths");}
				}
				prior.add(lo);
				prior.add(hi);
				seen++;
			}
		}
		check(seen==expected && cores.size()==expected && contexts.size()==expected && schedule.size()==expected,
			"All scheduled IDs must become exactly one core and context");
		try(Input input=new Input(path(fixture.root, prefix+".window_metadata.tsv"))){
			input.reader.nextLine();
			int rows=0;
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				final String id=p.parseString(0), source=p.parseString(1);
				final long genomeGc=source.equals("g8") ? 1990 : source.equals("g9") ? 796 : 0;
				check(p.parseLong(4)==1990 && p.parseLong(5)==genomeGc,
					"Exact GC counts must exclude five ambiguous bases on each of two contigs");
				int gc=0;
				for(byte base:cores.get(id)){if(base=='C' || base=='G'){gc++;}}
				check(p.parseInt(6)==gc && p.parseInt(7)==cores.get(id).length, "Window GC must describe the emitted core");
				rows++;
			}
			check(rows==expected, "Every emitted core needs its GC metadata");
		}
		int total=0, contributors=0, selected=0, min=Integer.MAX_VALUE, max=0;
		try(Input input=new Input(path(fixture.root, prefix+".strata.tsv"))){
			input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				p.set(row);
				check(p.parseInt(3)==p.parseInt(4), "Every occupied stratum must realize its fixed quota");
				selected+=p.parseInt(1);
				contributors+=p.parseInt(2);
				total+=p.parseInt(4);
				min=Math.min(min, p.parseInt(4));
				max=Math.max(max, p.parseInt(4));
			}
		}
		check(total==expected && max-min<=1 && contributors<=selected, "Balanced occupied-stratum quotas must conserve the target");
		if(expected==10){check(contributors<selected, "Zero-quota panel sources stay distinct from actual contributors");}
	}

	/** Owns small artificial source files; every intentional replacement stays in its fresh root. */
	private static final class Fixture{
		Fixture(final File directory){
			root=directory;
			genomes=path(root, "genomes.tsv");
			groups=path(root, "groups.tsv");
			ids=path(root, "ids.txt");
			loci=path(root, "loci.tsv");
			lengths=path(root, "lengths.fa");
		}
		void write(final boolean reverse){
			final ByteBuilder roster=new ByteBuilder(NcrnaCorpusSample.GENOME_HEADER+'\n');
			final ByteBuilder group=new ByteBuilder("genome\tassembly\ttaxid\tfasta_sha80\tsource_group\tpilot_exposed_group\n");
			final ByteBuilder mask=new ByteBuilder(NcrnaCmLocusAudit.HEADER+'\n'), panel=new ByteBuilder();
			for(int i=0; i<11; i++){
				final String name="g"+i, phylum=i==8 ? "Q" : i==9 ? "R" : "P";
				final String assembly=i==10 ? "GCF_0.2" : "GCF_"+i+".1";
				final int taxid=i==10 ? 1 : i+1;
				final byte[] sequence=new byte[1000];
				final String pattern=i==8 ? "C" : i==9 ? "GCAAA" : "A";
				for(int j=0; j<sequence.length; j++){sequence[j]=(byte)pattern.charAt(j%pattern.length());}
				Arrays.fill(sequence, 200, 205, (byte)'N');
				bases.put(name, sequence);
				roster.append(i).tab().append(name).tab().append(assembly).tab().append(taxid).append("\tbacteria\t").append(phylum)
					.append("\tunused\t").append(path(root, name+".fa")).append("\t1\t2000\t2\n");
				group.append(name).tab().append(assembly).tab().append(taxid).append("\tNA\tgroup").append(i==10 ? 0 : i).append("\tfalse\n");
				mask.append("RF00059_").append(i).tab().append(name).tab().append(assembly).tab().append(taxid).append("\tbacteria\t").append(phylum)
					.append("\ta\t+\t50-70\t21\t1\tx\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t0\t0\t1\n");
				if(i<10){panel.append(name).nl();}
			}
			NcrnaBackgroundSamplerTest.write(genomes, roster.toString());
			NcrnaBackgroundSamplerTest.write(groups, group.toString());
			NcrnaBackgroundSamplerTest.write(ids, panel.toString());
			NcrnaBackgroundSamplerTest.write(loci, mask.toString());
			writeLengths(lengths, new int[]{17, 18, 20, 23});
			writeSequences(reverse);
		}
		void writeSequences(final boolean reverse){
			for(int i=0; i<11; i++){writeGenome(path(root, "g"+i+".fa"), bases.get("g"+i), reverse);}
		}
		final File root;
		final String genomes, groups, ids, loci, lengths;
		final ObjectMap<String, byte[]> bases=new ObjectMap<String, byte[]>(String.class, byte[].class);
	}

	private static void writeGenome(final String path, final byte[] sequence, final boolean reverse){
		final ByteBuilder b=new ByteBuilder();
		b.append(reverse ? ">z\n" : ">a\n").append(sequence).nl();
		b.append(reverse ? ">a\n" : ">z\n").append(sequence).nl();
		write(path, b.toString());
	}

	private static void writeLengths(final String path, final int[] lengths){
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<lengths.length; i++){
			b.append('>').append(i).nl();
			for(int j=0; j<lengths[i]; j++){b.append('A');}
			b.nl();
		}
		write(path, b.toString());
	}

	/** Reverses records, retaining a schema header when present. */
	private static void reverse(final String source, final String target, final boolean header){
		final ArrayList<byte[]> rows=new ArrayList<byte[]>();
		try(Input input=new Input(source); Output output=new Output(target)){
			if(header){output.writer.println(input.reader.nextLine());}
			for(byte[] row; (row=input.reader.nextLine())!=null; ){rows.add(row);}
			Collections.reverse(rows);
			for(byte[] row:rows){output.writer.println(row); }
		}
	}

	/** Removes g10, the last source, leaving its g0 alias group deliberately incomplete. */
	private static void copyWithoutLast(final String source, final String target){
		try(Input input=new Input(source); Output output=new Output(target)){
			byte[] previous=input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null; ){
				output.writer.println(previous);
				previous=row;
			}
		}
	}

	private static void equalLines(final String left, final String right){
		try(Input a=new Input(left); Input b=new Input(right)){
			for(byte[] x=a.reader.nextLine(), y=b.reader.nextLine(); x!=null || y!=null; x=a.reader.nextLine(), y=b.reader.nextLine()){
				check(Arrays.equals(x, y), "Reordering equivalent input records changed output: "+left+" versus "+right);
			}
		}
	}

	private static int countRows(final String path){
		int count=0;
		try(Input input=new Input(path)){
			check(input.reader.nextLine()!=null, "Table needs a header: "+path);
			while(input.reader.nextLine()!=null){count++;}
		}
		return count;
	}

	/** Test-only replacement stays inside the fresh fixture directory. */
	private static void write(final String path, final String text){
		try(Output output=new Output(fileIO.FileFormat.testOutput(path, fileIO.FileFormat.TEXT, null, true, true, false, false))){
			output.writer.print(text);
		}
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void reject(final Runnable action, final String message){
		try{action.run();}catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(message), "Wrong diagnostic: "+failure.getMessage());
			return;
		}
		throw new AssertionError("Expected rejection: "+message);
	}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
