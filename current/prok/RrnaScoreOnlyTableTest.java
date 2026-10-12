package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/** Missing sparse keys can have a nonzero default when observed denominators differ. */
public final class RrnaScoreOnlyTableTest {
	public static void main(String[] args)throws Exception{
		final RrnaPositionalKmerTable.Table built=new RrnaPositionalKmerTable.Table("fixture","5prime",6,0,1.0/4096,2);
		final byte[] bases="ACGTACGTACGTACGTACGTACGTACGTACGTACGT".getBytes(StandardCharsets.US_ASCII);
		built.add(bases,10);built.add(bases,0);built.finish();
		final Path path=Files.createTempDirectory("score-only-table-").resolve("table.tsv.gz");built.write(path.toString());
		final RrnaPositionalKmerTable.Table full=RrnaPositionalKmerTable.load(path.toString(),"fixture"),compact=RrnaPositionalKmerTable.loadScores(path.toString(),"fixture");
		assert(compact.counts==null && full.counts!=null):"Inference must not allocate positional training counts";
		assert(full.scoreKey(0)!=0):"Fixture must exercise an omitted key with nonzero smoothed default";
		for(int word=0;word<full.space;word++){assert(Double.doubleToLongBits(full.scoreKey(word))==Double.doubleToLongBits(compact.scoreKey(word))):"Every observed and omitted key must retain exact double bits, key="+word;}
		final Path packedPath=path.resolveSibling("scores.tsv.gz");compact.writeScores(packedPath.toString());
		final RrnaPositionalKmerTable.Table packed=RrnaPositionalKmerTable.loadScores(packedPath.toString(),"fixture");
		assert(packed.counts==null && packed.records==full.records && packed.offset==full.offset && packed.radius==full.radius)
			: "Compact disk format must retain geometry and smoothing metadata without positional counts";
		for(int word=0; word<full.space; word++){
			assert(Double.doubleToLongBits(full.scoreKey(word))==Double.doubleToLongBits(packed.scoreKey(word)))
				: "Compact shipping rows and omitted defaults must reproduce every original double bit, key="+word;
		}
		boolean rejected=false;try{RrnaPositionalKmerTable.load(packedPath.toString(),"fixture");}catch(IllegalArgumentException e){rejected=true;}
		assert(rejected) : "Scores-only resources cannot masquerade as full training-count tables";
		System.out.println("RRNA_SCORE_ONLY_TABLE_FIXTURE_PASS keys="+full.space);
	}
}
