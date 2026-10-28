package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CobalaminGenomicContext;
import prok.NcrnaReviewJoin;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Tests strand, circular topology, overlap units, ties and intact product text.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminGenomicContextTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);final String gff=root.resolve("source.gff").toString(), hash="00000000000000000000";
		final String[] ids={"plus", "minus", "circular", "linear"};final int[] starts={100,800,980,980}, stops={119,819,999,999}, counts={2,1,0,0};
		final StringBuilder loci=new StringBuilder(NcrnaReviewUnion.LOCUS_HEADER).append('\n'), context=new StringBuilder(NcrnaReviewJoin.CONTEXT_HEADER).append('\n');
		for(int i=0; i<ids.length; i++){
			final String strand=i==1?"-":"+";
			loci.append(ids[i]).append("\tg\tA\tgenome.fa\t").append(ids[i]).append("\tcobalamin\t").append(strand).append('\t').append(starts[i]).append('\t').append(stops[i])
				.append("\t1000\t0\t0\t20\tPREDICTION_UNRESOLVED\n");
			context.append(ids[i]).append("\tg\tA\t").append(ids[i]).append('\t').append(starts[i]).append('\t').append(stops[i]).append('\t').append(strand)
				.append("\t20\t189\t-169\t0.105820\t").append(counts[i]).append("\t1\t0\t0\t").append(i==0?1:0).append('\t').append(i<2?1:0)
				.append("\t0\t").append(i<2?1:0).append('\t').append(i==0?1:0).append("\t0\t").append(counts[i]>0?"FEATURE_OVERLAP":"INTERGENIC")
				.append('\t').append(gff).append('\t').append(hash).append("\tgenome.fa\t").append(hash).append('\n');
		}
		write(root.resolve("loci.tsv"), loci.toString());write(root.resolve("context.tsv"), context.toString());
		final String source=feature("plus", "CDS", 105, 110, "-", "gene=x;product=opposite overlap;")+
			feature("plus", "tRNA", 110, 130, "+", "product=tRNA-Arg;")+
			feature("plus", "CDS", 150, 170, "-", "gene=x;product=closer opposite gene;")+
			feature("plus", "CDS", 200, 230, "+", "gene=x;product=other enzyme with spaces;")+
			feature("plus", "CDS", 200, 249, "+", "gene=metH;product=Cobalamin dependent methionine synthase;")+
			feature("minus", "CDS", 790, 830, "-", "gene=x;product=containing CDS;")+
			feature("minus", "CDS", 700, 749, "-", "gene=btuF;product=vitamin B12 transport protein;")+
			feature("minus", "CDS", 750, 780, "+", "gene=x;product=closer opposite gene;")+
			feature("circular", "CDS", 10, 50, "+", "gene=cobA;product=precorrin methyltransferase;")+
			feature("linear", "CDS", 10, 50, "+", "gene=cobA;product=precorrin methyltransferase;")+
			// Topology deliberately follows CDS rows; source ordering must not change selection.
			feature("plus", "region", 1, 1000, "+", "Is_circular=false;")+feature("minus", "region", 1, 1000, "+", "Is_circular=false;")+
			feature("circular", "region", 1, 1000, "+", "Is_circular=true;")+feature("linear", "region", 1, 1000, "+", "Is_circular=false;");
		write(root.resolve("source.gff"), source);
		CobalaminGenomicContext.main(new String[]{"loci="+root.resolve("loci.tsv"), "context="+root.resolve("context.tsv"), "outprefix="+root.resolve("result"), "expected=4"});
		final ByteFile in=RrnaResourceIO.open(root.resolve("result.neighborhood.tsv").toString());final LineParser1 p=new LineParser1('\t');int rows=0;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);rows++;check(p.terms()==18, "Neighborhood schema must retain both nearest same-strand and any-strand CDS");
			if(p.termEquals("plus", 0)){
				check(p.termEquals("CDS_OPPOSITE+tRNA", 1) && p.parseInt(4)==6 && p.parseDouble(5)==.3 && p.termEquals("false", 6), "CDS overlap fraction uses call length and retains RNA co-overlap");
				check(p.parseInt(7)==80 && p.termEquals("true", 8) && p.parseInt(11)==2 && p.termEquals("other enzyme with spaces", 10), "Tied nearest genes retain full raw product and any-tied B12 hint");
				check(p.parseInt(12)==30 && p.termEquals("-", 17), "Closer opposite-strand CDS must not replace the same-strand candidate");
			}else if(p.termEquals("minus", 0)){
				check(p.termEquals("CDS_SAME", 1) && p.parseInt(2)==20 && p.termEquals("true", 6) && p.parseInt(7)==50 && p.parseInt(12)==19,
					"Minus-strand direction reverses neighbor search; a containing CDS remains overlap, not downstream");
			}else if(p.termEquals("circular", 0)){check(p.termEquals("INTERGENIC", 1) && p.parseInt(7)==10, "Explicit circular topology permits an origin-crossing downstream neighbor");}
			else{check(p.termEquals("linear", 0) && p.termEquals("NA", 7), "The same coordinate configuration must not wrap on a linear contig");}
		}}finally{check(!in.close(), "Fixture neighborhood read failed");}
		check(rows==4, "Every frozen prediction must appear once");
		check(CobalaminGenomicContext.downstreamGap(100,119,'+',120,150,1000,false)==0
			&& CobalaminGenomicContext.downstreamGap(800,819,'-',770,799,1000,false)==0, "Adjacent CDS gap is zero on both strands");
		check(CobalaminGenomicContext.downstreamGap(1,20,'-',950,980,1000,true)==20, "Minus-strand circular distance crosses the origin once");
		check(CobalaminGenomicContext.hint("methylase", "generic enzyme").equals("NONE") && CobalaminGenomicContext.hint("cbiA", "enzyme").equals("gene:cbi"), "Specific metH name must not match arbitrary methylase text; cbi family prefix is intentional");
		check(CobalaminGenomicContext.attribute("product=full product with spaces;gene=x;", "product").equals("full product with spaces"), "GFF product values end at semicolons, not whitespace");
		System.out.println("COBALAMIN_GENOMIC_CONTEXT_TEST_PASS");
	}
	private static String feature(String id, String type, int a, int z, String strand, String attrs){return id+"\tsource\t"+type+'\t'+a+'\t'+z+"\t.\t"+strand+"\t.\t"+attrs+'\n';}
	private static void write(Path path, String text){final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
