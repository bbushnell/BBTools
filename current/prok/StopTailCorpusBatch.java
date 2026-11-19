package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Runs independent annotated genomes serially within one JVM, preserving native extraction.
 * Separate processes may handle separate manifests; GFF parser globals are never shared concurrently.
 * Defaults retain whole SHORT proteins; explicit minlen/cap/maxoriginal also support LONG extraction.
 * @author Fischl
 */
public final class StopTailCorpusBatch{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, StopTailCorpusBatch.class, false);
		try{
			final HashMap<String,String> options=new HashMap<String,String>();
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				require(pair.length==2 && !pair[1].isEmpty() && (pair[0].equals("manifest") || pair[0].equals("manifestsha80")
					|| pair[0].equals("out") || pair[0].equals("rejectpseudo")
					|| pair[0].equals("reject5primepartial") || pair[0].equals("minlen")
					|| pair[0].equals("cap") || pair[0].equals("maxoriginal")),
					"Batch expects manifest, out and optional rejection/length policies; legacy manifestsha80 is ignored");
				require(options.put(pair[0],pair[1])==null, "Duplicate batch option");
			}
			require(options.containsKey("manifest") && options.containsKey("out"),
				"Batch needs explicit manifest and out; input identity is the committed literal path list");
			run(Paths.get(options.get("manifest")), Paths.get(options.get("out")),
				Parse.parseBoolean(options.getOrDefault("rejectpseudo", "f")),
				Parse.parseBoolean(options.getOrDefault("reject5primepartial", "f")),
				Integer.parseInt(options.getOrDefault("minlen", "40")),
				Integer.parseInt(options.getOrDefault("cap", "200")),
				Integer.parseInt(options.getOrDefault("maxoriginal", "99")));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Calls the existing extractor once per listed genome; legacy pin columns are metadata only. */
	private static void run(final Path manifest, final Path out,
			final boolean rejectPseudo, final boolean rejectFivePrime,
			final int minlen, final int cap, final int maxOriginal) throws Exception{
		require(minlen>=1 && cap>=minlen && maxOriginal>=minlen,
			"Require 1<=minlen<=cap and maxoriginal>=minlen; use2147483647 for no upper length limit");
		final String label=minlen==40 && cap==200 && maxOriginal==99 ? "SHORT_CORPUS" : "STOP_TAIL_CORPUS";
		requireInput(manifest);
		require(!Files.exists(out), "Batch output must be fresh; failed batches are preserved");
		Shared.setThreads(1);
		ReadWrite.USE_UNPIGZ=false;
		Files.createDirectory(out);
		Files.copy(manifest, out.resolve("inputs.tsv"));
		final HashSet<String> seen=new HashSet<String>();
		final long[] totals=new long[METRICS.length];
		final ByteFile input=ByteFile.makeByteFile1(manifest.toString(), false);
		final LineParser1 row=new LineParser1('\t');
		try{
			final byte[] header=input.nextLine();
			final boolean legacy=java.util.Arrays.equals(header, LEGACY_HEADER.getBytes(StandardCharsets.US_ASCII));
			require(legacy || java.util.Arrays.equals(header, HEADER.getBytes(StandardCharsets.US_ASCII)),
				"Expected path-only or legacy extraction manifest header");
			final int offset=legacy ? 2 : 0;
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				row.set(line);
				require(row.terms()==10+offset, "Extraction manifest row differs from its declared schema");
				final String genome=row.parseString(0), source=row.parseString(1);
				require(genome.matches("[A-Za-z0-9][A-Za-z0-9_.-]*") && seen.add(genome), "Unique path-safe genome IDs required");
				final Path fasta=Paths.get(row.parseString(2)), gff=Paths.get(row.parseString(legacy ? 4 : 3));
				requireInput(fasta); requireInput(gff);
				final Path genomeOut=out.resolve(genome);
				Files.createDirectory(genomeOut);
				StopTailCorpus.main(new String[]{"in="+fasta,"gff="+gff,"genome="+genome,"source="+source,
					"taxid="+row.parseString(4+offset), "phylumtaxid="+row.parseString(5+offset),
					"taxnamespace="+row.parseString(6+offset),
					"split="+row.parseString(7+offset),"table="+row.parseString(8+offset),"tableassumed="+row.parseString(9+offset),
					"out="+genomeOut.resolve("tails.fa.gz"),"outmeta="+genomeOut.resolve("tails.tsv.gz"),
					"outstats="+genomeOut.resolve("stats.tsv"),"outaudit="+genomeOut.resolve("audit.tsv.gz"),
					"rejectpseudo="+rejectPseudo,"reject5primepartial="+rejectFivePrime,
					"minlen="+minlen,"cap="+cap,"maxoriginal="+maxOriginal,"overwrite=f","t=1","zl=9","bgzip=t"});
				addStats(genomeOut.resolve("stats.tsv"), genome, totals, maxOriginal==Integer.MAX_VALUE);
				System.err.println(label+"_PROGRESS genomes="+seen.size()+" emitted="+totals[1]);
			}
		}finally{require(!input.close(), "Batch manifest read failed");}
		require(!seen.isEmpty(), "An empty batch cannot establish extraction completion");
		final ByteBuilder summary=new ByteBuilder("metric\tcount\ngenomes\t").append(seen.size()).nl();
		for(int i=0; i<METRICS.length; i++){summary.append(METRICS[i]).tab().append(totals[i]).nl();}
		final ByteStreamWriter writer=new ByteStreamWriter(out.resolve("summary.tsv").toString(), false, false, false);
		writer.start(); writer.print(summary);
		require(!writer.poisonAndWait(), "Batch summary write failed");
		Files.write(out.resolve("PASS"), new ByteBuilder(label).append("_BATCH_PASS genomes=").append(seen.size())
			.append(" emitted=").append(totals[1]).nl().toBytes());
		System.out.print(summary.toString());
	}

	/** Sum native counters once; do not infer missing values from output file size. */
	private static void addStats(final Path path, final String genome, final long[] totals, final boolean unlimited){
		final boolean[] seen=new boolean[METRICS.length];
		final ByteFile input=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 row=new LineParser1('\t');
		try{
			require(java.util.Arrays.equals(input.nextLine(),
				"schema\tgenome_id\tmetric\tcount".getBytes(StandardCharsets.US_ASCII)),
				"Expected native extraction counters");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				row.set(line);
				require(row.terms()==4 && row.termEquals("stop_tail_v1",0) && row.termEquals(genome,1),
					"Counter provenance differs");
				final long count=row.parseLong(3);
				require(count>=0, "Extraction counts must be nonnegative");
				for(int i=0; i<METRICS.length; i++){
					if(row.termEquals(METRICS[i],2)){
						require(!seen[i], "Each native counter must occur once"); seen[i]=true;
						require(!unlimited || !METRICS[i].equals("too_long") || count==0,
							"An unlimited original-length range cannot reject a gene as too long");
						totals[i]=Math.addExact(totals[i],count);
					}
				}
			}
		}finally{require(!input.close(), "Native counter read failed");}
		for(int i=0; i<METRICS.length; i++){
			// StopTailCorpus.write omits this zero counter only when no upper bound is applied.
			require(seen[i] || unlimited && METRICS[i].equals("too_long"), "A required native counter is missing: "+METRICS[i]);
		}
	}
	/** Metadata-only path validation; input bytes are not hashed or checked before/after. */
	private static void requireInput(final Path path){
		require(path.isAbsolute() && Files.isRegularFile(path), "An explicit regular input path is required: "+path);
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
	public static final String HEADER="genome_id\tsource\tfasta\tgff\ttaxid\tphylum_taxid"
		+"\ttaxonomy_namespace\tsplit\ttable\ttable_assumed";
	/** Accepted only for compatibility; both SHA columns are ignored, never verified. */
	public static final String LEGACY_HEADER="genome_id\tsource\tfasta\tfasta_sha80\tgff\tgff_sha80"
		+"\ttaxid\tphylum_taxid\ttaxonomy_namespace\tsplit\ttable\ttable_assumed";
	private static final String[] METRICS={"cds_total","emitted","duplicate_alias","too_short","too_long",
		"internal_stop","ambiguous_residue","tail_residues",
		"pseudo","five_prime_partial","three_prime_partial"};
}
