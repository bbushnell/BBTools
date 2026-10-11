package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import fileIO.ReadWrite;
import shared.Shared;

/** Exact text/count/probability parity and malformed-input tests for nucleotide HBMC.
 * @author Keqing
 */
public final class NucleotideHbmCodecTest {
	public static void main(String[] args) throws Exception{
		Shared.setThreads(1);
		if(args.length==2 && args[0].startsWith("in=") && args[1].startsWith("other=")){
			compare(TrnaConsensusBuilder.loadModels(args[0].substring(3)),TrnaConsensusBuilder.loadModels(args[1].substring(6)));
			return;
		}
		if(args.length!=1 || !args[0].startsWith("dir=")){throw new IllegalArgumentException("Use dir=fresh_directory or in=dense.hbm other=compact.hbmc");}
		final Path dir=java.nio.file.Paths.get(args[0].substring(4));Files.createDirectory(dir);
		final String source="#HBM\t1\n>first n=3\tdepth=999\naCgTNACG\n"
			+"3\t0\t0\t0\t1\n0\t64\t1\t0\t0\n0\t0\t4096\t0\t0\n0\t0\t0\t63\t0\n"
			+"0\t0\t0\t0\t0\n1000000\t0\t1\t0\t0\n1\t2\t3\t4\t5\n0\t0\t8\t0\t2\n"
			+">second\tdepth=12\nACGTACGT\n12\t0\t0\t0\t0\n0\t12\t0\t0\t0\n0\t0\t12\t0\t0\n0\t0\t0\t12\t0\n"
			+"12\t0\t0\t0\t0\n0\t12\t0\t0\t0\n0\t0\t12\t0\t0\n0\t0\t0\t12\t0\n";
		final Path dense=write(dir,"input.hbm",source), compact=dir.resolve("output.hbmc"), restored=dir.resolve("restored.hbm");
		NucleotideHbmCodec.convert(dense.toString(),compact.toString());
		NucleotideHbmCodec.convert(compact.toString(),restored.toString());
		check(Arrays.equals(Files.readAllBytes(dense),Files.readAllBytes(restored)),"Headers, sequence case, counts and order must round-trip exactly");
		compare(TrnaConsensusBuilder.loadModels(dense.toString()),TrnaConsensusBuilder.loadModels(compact.toString()));
		final Path compressed=dir.resolve("output.hbmc.gz"), decoded=dir.resolve("decoded.hbm.gz");
		NucleotideHbmCodec.convert(dense.toString(),compressed.toString());
		NucleotideHbmCodec.convert(compressed.toString(),decoded.toString());
		check(ReadWrite.readString(decoded.toString()).equals(source),"Compressed conversion preserves decoded text");
		compare(TrnaConsensusBuilder.loadModels(dense.toString()),TrnaConsensusBuilder.loadModels(compressed.toString()));
		final String oldDori=dna.Data.DORI_RESOURCES, oldNersc=dna.Data.NERSC_RESOURCES;
		try{
			dna.Data.DORI_RESOURCES=dir.toString()+"/";
			final String probe="nucleotide_codec_resource_fixture.hbm";
			check(CallGenes.defaultNcrnaModel(probe).equals(probe),"Missing compact model retains legacy resource name");
			final Path packed=dir.resolve(probe+"c.gz");Files.copy(compressed,packed);
			check(CallGenes.defaultNcrnaModel(probe).equals(packed.toString()),"Compact default resolves gzip resource");
			final Path plain=dir.resolve(probe+"c");Files.copy(compact,plain);
			check(CallGenes.defaultNcrnaModel(probe).equals(plain.toString()),"Exact compact name precedes its gzip alternative");
			final Path later=Files.createDirectory(dir.resolve("later"));dna.Data.NERSC_RESOURCES=later.toString()+"/";
			Files.copy(compressed,later.resolve("codec_priority_fixture.hbmc.gz"));
			Files.copy(dense,dir.resolve("codec_priority_fixture.hbm"));
			check(CallGenes.defaultNcrnaModel("codec_priority_fixture.hbm").equals("codec_priority_fixture.hbm"),"Later compact location must not shadow an earlier legacy model");
			final Path cwd=Files.createTempFile(java.nio.file.Paths.get("."),"hbmc-cwd-",".hbm");
			try{
				final String name=cwd.getFileName().toString();final Path alternate=dir.resolve(name+"c.gz");Files.copy(compressed,alternate);
				check(CallGenes.defaultNcrnaModel(name).equals(name),"Generic ncRNA keeps an existing CWD model");
				check(CallGenes.defaultNcrnaModel(name,false).equals(alternate.toString()),"tRNA retains its historical resource-root-only lookup");
			}finally{Files.delete(cwd);}
		}finally{dna.Data.DORI_RESOURCES=oldDori;dna.Data.NERSC_RESOURCES=oldNersc;}
		final String valid=new String(Files.readAllBytes(compact),StandardCharsets.UTF_8);
		reject(dir,"no-footer.hbmc",valid.substring(0,valid.indexOf("#end\t")));
		reject(dir,"bad-footer.hbmc",valid.replace("#end\t2\t16","#end\t1\t16"));
		reject(dir,"trailing.hbmc",valid+"junk\n");
		reject(dir,"short.hbmc","#HBMC_NT\t1\n>x\tdepth=1\nAC\n1\n#end\t1\t2\n");
		reject(dir,"empty.hbmc","#HBMC_NT\t1\n>x\nA\n\n#end\t1\t1\n");
		reject(dir,"bad-digit.hbmc","#HBMC_NT\t1\n>x\nA\np\n#end\t1\t1\n");
		reject(dir,"overflow.hbmc","#HBMC_NT\t1\n>x\nA\noooooo\n#end\t1\t1\n");
		rejectDense(dir,"negative.hbm","#HBM\t1\n>x\nA\n-1\t0\t0\t0\t0\n");
		rejectDense(dir,"depth-overflow.hbm","#HBM\t1\n>x\nA\n2147483647\t1\t0\t0\t0\n");
		final Path input=write(dir,"builder.fa",">tRNA-Test(AAA)\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
		BaseGraph[] reference=null;
		for(String extension:new String[]{"hbm","hbm.gz","hbmc","hbmc.gz"}){
			final Path model=dir.resolve("builder."+extension);
			TrnaConsensusBuilder.main(new String[]{"in="+input,"out="+dir.resolve("consensus."+extension+".fa"),
				"outmodel="+model,"cluster=f","minclustersize=1","mingroupsize=1","t=1"});
			final BaseGraph[] found=TrnaConsensusBuilder.loadModels(model.toString());
			if(reference==null){reference=found;}else{compare(reference,found);}
		}
		for(String extension:new String[]{"hbm","hbmc"}){
			final Path model=dir.resolve("builder."+extension);
			TrnaConsensusBuilder.main(new String[]{"in="+input,"out="+dir.resolve("overwrite."+extension+".fa"),
				"outmodel="+model,"cluster=f","minclustersize=1","mingroupsize=1","t=1","ow=t"});
			compare(reference,TrnaConsensusBuilder.loadModels(model.toString()));
			final byte[] before=Files.readAllBytes(model);boolean refused=false;
			try{TrnaConsensusBuilder.main(new String[]{"in="+input,"out="+dir.resolve("refused."+extension+".fa"),
				"outmodel="+model,"cluster=f","minclustersize=1","mingroupsize=1","t=1","ow=f"});}
			catch(RuntimeException | AssertionError expected){refused=true;}
			check(refused && Arrays.equals(before,Files.readAllBytes(model)),"ow=f must refuse without modifying "+extension);
		}
		System.out.println("NUCLEOTIDE_HBMC_SELFTEST_PASS malformed_cases=9 builder_formats=4 overwrite_formats=2 cwd_priority=exact");
	}

	/** Every reconstructed count/probability is checked; alignment probes cover ten models. */
	static void compare(BaseGraph[] a,BaseGraph[] b){
		check(a.length==b.length,"Model count");long positions=0;int probes=0;
		for(int m=0;m<a.length;m++){
			final BaseGraph x=a[m],y=b[m];
			check(x.name.equals(y.name) && x.numericID==y.numericID && Arrays.equals(x.original,y.original),"Model identity/order at "+m);
			check(x.ref.length==y.ref.length,"Reference length at "+m);
			for(int p=0;p<x.ref.length;p++){
				check(Arrays.equals(BaseGraphHelper.getCounts(x,p),BaseGraphHelper.getCounts(y,p)),"Counts at "+m+":"+p);
				check(Arrays.equals(x.ref[p].acgtWeight,y.ref[p].acgtWeight) && Arrays.equals(x.ref[p].acgtProb,y.ref[p].acgtProb),"Weights/probabilities at "+m+":"+p);
				check(x.ref[p].countSum==y.ref[p].countSum && x.ref[p].weightSum==y.ref[p].weightSum &&
					x.del[p].countSum==y.del[p].countSum && x.del[p].weightSum==y.del[p].weightSum,"Depths at "+m+":"+p);
				positions++;
			}
			if(m<10){
				final byte[] mutated=x.original.clone();mutated[mutated.length/2]=mutated[mutated.length/2]=='A'?(byte)'C':(byte)'A';
				for(byte[] query:new byte[][]{x.original,mutated}){
					final float sx=TrnaConsensusBuilder.scoreAgainstModel(query,x),sy=TrnaConsensusBuilder.scoreAgainstModel(query,y);
					check(Float.floatToIntBits(sx)==Float.floatToIntBits(sy),"Alignment score at model "+m);probes++;
				}
			}
		}
		System.out.println("NUCLEOTIDE_HBM_PARITY_PASS models="+a.length+" positions="+positions+" score_probes="+probes);
	}
	private static Path write(Path dir,String name,String text) throws Exception{
		final Path path=dir.resolve(name);Files.write(path,text.getBytes(StandardCharsets.UTF_8));return path;
	}
	private static void reject(Path dir,String name,String text) throws Exception{
		final Path path=write(dir,name,text);boolean failed=false;
		try{NucleotideHbmCodec.loadModels(path.toString());}catch(IllegalArgumentException expected){failed=true;}
		check(failed,"Malformed compact input accepted: "+name);
	}
	private static void rejectDense(Path dir,String name,String text) throws Exception{
		final Path path=write(dir,name,text);boolean failed=false;
		try{NucleotideHbmCodec.convert(path.toString(),dir.resolve(name+".hbmc").toString());}catch(IllegalArgumentException expected){failed=true;}
		check(failed,"Malformed dense input accepted: "+name);
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	private NucleotideHbmCodecTest(){}
}
