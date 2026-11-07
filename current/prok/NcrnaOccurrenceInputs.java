package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import map.IntHashMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaObservationData.Symbols;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/** Prepares source-local extraction metadata for a complete occurrence population.
 * Forked from template.A_SampleByteFile. No sequence reader, admission filter or
 * biological annotation label is introduced. One metadata thread owns all state.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOccurrenceInputs{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaOccurrenceInputs tool=new NcrnaOccurrenceInputs(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	/** Parses native options; all outputs belong to a new directory. */
	public NcrnaOccurrenceInputs(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		try{
			final Parser parser=new Parser();
			for(final String arg:args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}
				if(a.equals("roster")){roster=b;}
				else if(a.equals("flank")){flank=Parse.parseIntKMG(b);}
				else if(a.equals("expected")){expected=Parse.parseIntKMG(b);}
				else if(parser.parse(arg, a, b)){
					// Shared BBTools options.
				}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			in=parser.in1;
			out=parser.out1;
			if(in==null || roster==null || out==null || flank<0 || expected<1 || new File(out).exists()){
				throw new IllegalArgumentException("Required: in=occurrences.tsv roster=genomes.tsv out=fresh_directory expected=N flank>=0");
			}
			if(parser.append || parser.maxReads>=0 || parser.samplerate!=1 || parser.in2!=null || parser.out2!=null){
				throw new IllegalArgumentException("Occurrence extraction metadata requires the complete unpaired population");
			}
		}catch(RuntimeException | Error failure){
			Shared.closeStream(outstream);
			throw failure;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Input preparation      ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		assert(sources.isEmpty() && features.isEmpty()) : "One preparer binds one admitted population";
		readRoster();
		readOccurrences();
		if(!new File(out).mkdir() || !new File(out, "gff").mkdir()){
			throw new IllegalArgumentException("Cannot create fresh extraction directory: "+out);
		}
		writeSources();
		writeFeatures();
		timer.stop();
		outstream.println("NCRNA_OCCURRENCE_INPUTS_PASS occurrences="+features.size()+" sources="+positiveSources
			+" expected_context_bases="+contextBases+" time="+timer);
	}

	/** Roster identity includes the exact assembly and paired FASTA, not TaxID alone. */
	private void readRoster(){
		assert(sources.isEmpty()) : "Source indices must be established before occurrence joins";
		try(Input input=new Input(roster)){
			header(input, NcrnaLearningPopulation.ROSTER_HEADER);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				parse(row, 11);
				final int index=parser.parseInt(0), genome=symbol(1), assembly=symbol(2), taxid=parser.parseInt(3);
				final int domain=symbol(4), phylum=symbol(5), fasta=symbol(7);
				require(index==sources.size() && byGenome.get(genome)<0 && taxid>0, "Roster indices/names/TaxIDs must be unique, dense and valid");
				identifier(genome);
				final String path=parser.parseString(7);
				require(new File(path).isAbsolute() && path.indexOf('=')<0,
					"CutGff requires an absolute FASTA path without '='; it checks file access when extraction starts");
				byGenome.put(genome, index);
				sources.add(new Source(index, genome, assembly, taxid, domain, phylum, fasta, row));
			}
		}
		require(!sources.isEmpty(), "A complete source roster is required");
	}

	/** Calculates only exact inclusive coordinate transforms; no occurrence is filtered. */
	private void readOccurrences(){
		assert(!sources.isEmpty()) : "Every admitted occurrence must join an established source";
		final IntHashMap seen=new IntHashMap();
		try(Input input=new Input(in)){
			header(input, OCCURRENCE_HEADER);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				parse(row, 12);
				final int id=symbol(0), sourceIndex=byGenome.get(symbol(1));
				require(sourceIndex>=0 && seen.get(id)<0, "Occurrence IDs must be unique and every source must be in the roster");
				identifier(id);
				seen.put(id, features.size());
				final Source source=sources.get(sourceIndex);
				require(symbol(2)==source.assembly && parser.parseInt(3)==source.taxid, "Occurrence assembly/TaxID differs from its exact source");
				final int seqid=symbol(4), start=parser.parseInt(6), stop=parser.parseInt(7), length=parser.parseInt(11);
				identifier(seqid);
				require(start>=1 && stop>=start && stop<=length, "Occurrence bounds must fit the declared source contig");
				final boolean minus=parser.termEquals('-', 5);
				require(minus || parser.termEquals('+', 5), "An explicit genomic strand is required");
				final int left=Math.min(flank, start-1), right=Math.min(flank, length-stop);
				final Feature feature=new Feature(id, source, seqid, start, stop, length, minus, minus ? right : left, minus ? left : right);
				features.add(feature);
				source.features.add(feature);
				contextBases+=(long)stop-start+1+left+right;
			}
		}
		require(features.size()==expected, "Occurrence count differs from the explicitly declared full population");
	}

	/** Emits one GFF per positive source and an ordered paired-input configuration. */
	private void writeSources(){
		assert(features.size()==expected) : "Output creation follows complete population validation";
		try(Output config=output("sources.config"); Output rosterOutput=output("genomes.tsv")){
			rosterOutput.writer.println(NcrnaLearningPopulation.ROSTER_HEADER);
			for(final Source source:sources){
				if(!source.features.isEmpty()){
					positiveSources++;
					final String gff=new File(out, "gff/"+source.index+".gff").getAbsolutePath();
					config.writer.print(buffer.clear().append("in=").append(symbols.bytes(source.fasta)).nl()
						.append("gff=").append(gff).nl());
					try(Output gffOutput=new Output(gff)){
						gffOutput.writer.println("##gff-version 3");
						for(final Feature feature:source.features){
							buffer.clear().append(symbols.bytes(feature.seqid)).append("\toccurrence\tncRNA\t")
								.append(feature.start).tab().append(feature.stop).append("\t.\t")
								.append(feature.minus ? '-' : '+').append("\t.\t");
							attributes(feature);
							gffOutput.writer.print(buffer.nl());
						}
					}
				}
				parse(source.row, 11);
				buffer.clear();
				for(int i=0; i<11; i++){
					if(i>0){buffer.tab();}
					if(i==8){buffer.append(source.features.size());}
					else{
						parser.setBounds(i);
						buffer.append(source.row, parser.a(), parser.b()-parser.a());
					}
				}
				rosterOutput.writer.print(buffer.nl());
			}
		}
	}

	/** Separate record and occurrence IDs preserve both CutGff headers and admission joins. */
	private void writeFeatures(){
		assert(positiveSources>0) : "A nonempty validated population must have at least one extraction source";
		try(Output loci=output("occurrence_loci.tsv"); Output cores=output("core.gff");
				Output coverage=output("coverage_split.tsv"); Output summary=output("summary.tsv")){
			loci.writer.println(NcrnaLearningPopulation.OCCURRENCE_HEADER);
			cores.writer.println("##gff-version 3");
			coverage.writer.println("locus_id\tgenome\tsource_group\tstratum\trepresentative_id\tfold\tpilot_exposed_group");
			for(final Feature f:features){
				final Source s=f.source;
				buffer.clear().append("ID=").append(symbols.bytes(f.id)).tab().append(symbols.bytes(f.id)).tab()
					.append(symbols.bytes(s.genome)).tab().append(symbols.bytes(s.assembly)).tab().append(s.taxid)
					.tab().append(symbols.bytes(s.domain)).tab().append(symbols.bytes(s.phylum)).tab()
					.append(symbols.bytes(f.seqid)).tab().append(f.minus ? '-' : '+').tab().append(f.start).tab()
					.append(f.stop).tab().append(f.stop-f.start+1).tab().append(f.contigLength).tab().append(f.left).tab().append(f.right);
				loci.writer.print(buffer.nl());
				buffer.clear().append("ID=").append(symbols.bytes(f.id)).append("\toccurrence\tncRNA\t")
					.append(f.left+1).tab().append(f.left+f.stop-f.start+1).append("\t.\t+\t.\t");
				attributes(f);
				cores.writer.print(buffer.nl());
				buffer.clear().append("ID=").append(symbols.bytes(f.id)).tab().append(symbols.bytes(s.genome))
					.append("\ttid_").append(s.taxid).append("\tALL_ADMITTED_IN_CORPUS\tID=").append(symbols.bytes(f.id))
					.append("\tTRAIN\tNA\n");
				coverage.writer.print(buffer);
			}
			summary.writer.println("occurrences\tpositive_sources\troster_sources\tflank\texpected_context_bases");
			summary.writer.print(buffer.clear().append(features.size()).tab().append(positiveSources).tab()
				.append(sources.size()).tab().append(flank).tab().append(contextBases).nl());
		}
	}

	private void attributes(final Feature feature){
		buffer.append("ID=").append(symbols.bytes(feature.id)).append(" source=").append(symbols.bytes(feature.source.genome))
			.append(" tid_").append(feature.source.taxid);
	}
	private Output output(final String name){return new Output(new File(out, name).getPath());}
	private void parse(final byte[] row, final int width){
		current=row;
		parser.set(row);
		require(parser.terms()==width, "Unexpected metadata table width");
	}
	private int symbol(final int field){
		parser.setBounds(field);
		require(parser.b()>parser.a(), "Empty metadata identity");
		return symbols.intern(current, parser.a(), parser.b());
	}
	/** Tokens used in GFF/header identity fields cannot contain separators or controls. */
	private void identifier(final int id){
		for(final byte b:symbols.bytes(id)){
			require(b>32 && b<127 && b!=';' && b!='=' && b!='%', "Unsafe GFF/header identity token");
		}
	}
	private static void header(final Input input, final String expected){
		final byte[] row=input.reader.nextLine();
		require(row!=null && expected.equals(new String(row, java.nio.charset.StandardCharsets.UTF_8)), "Unexpected metadata table header");
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Retained data        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(int index_, int genome_, int assembly_, int taxid_, int domain_, int phylum_, int fasta_, byte[] row_){
			index=index_;
			genome=genome_;
			assembly=assembly_;
			taxid=taxid_;
			domain=domain_;
			phylum=phylum_;
			fasta=fasta_;
			row=row_;
		}
		final int index, genome, assembly, taxid, domain, phylum, fasta;
		final byte[] row;
		final ArrayList<Feature> features=new ArrayList<Feature>();
	}
	private static final class Feature{
		Feature(int id_, Source source_, int seqid_, int start_, int stop_, int contigLength_, boolean minus_, int left_, int right_){
			id=id_;
			source=source_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			contigLength=contigLength_;
			minus=minus_;
			left=left_;
			right=right_;
		}
		final int id, seqid, start, stop, contigLength, left, right;
		final Source source;
		final boolean minus;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, roster, out;
	private int flank=50, expected=-1, positiveSources;
	private long contextBases;
	private PrintStream outstream;
	private final Symbols symbols=new Symbols();
	private final IntHashMap byGenome=new IntHashMap();
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Feature> features=new ArrayList<Feature>();
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private byte[] current;
	private static final String OCCURRENCE_HEADER="occurrence_id\tgenome\tassembly\ttaxid\tseqid\tstrand\tstart1\tstop1\tmodel_id\tmodel_sha80\trepresentative_physical_id\tcontig_length";
}
