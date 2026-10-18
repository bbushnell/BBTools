package test.ml;

import ml.RegressionTrainer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.stream.Stream;

/**
 * Constructor-only RegressionTrainer validation checks for netin= without dims=.
 * Each invocation constructs one trainer and never calls process(). The unused
 * netin= path is deliberately temp-relative: construction only checks in= and
 * out= paths, while loadInitialNet() would load netin= later during process().
 * @author Jean
 */
public final class RegressionTrainerConstructorValidationTest{

	/** Runs one named constructor validation case. */
	public static void main(String[] args) throws Exception{
		boolean enabled=false; assert(enabled=true);
		if(!enabled || args.length!=1){throw new IllegalArgumentException("Run with -ea and one case name");}
		final String name=args[0];
		final Path dir=Files.createTempDirectory("regressiontrainer-constructor-validation-");
		try{
			final Path input=dir.resolve("rows.tsv");
			Files.write(input, ROWS.getBytes(StandardCharsets.US_ASCII));
			runCase(name, input, dir.resolve("out.bbnet"), dir.resolve("unused-net.bbnet"));
			System.out.println("RegressionTrainerConstructorValidationTest PASS "+name);
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** Dispatches the requested constructor-only case. */
	private static void runCase(String name, Path input, Path output, Path netIn){
		if(name.equals("valid-netin")){
			construct(input, output, netIn, "epochs=1", "batch=1");
		}else if(name.equals("missing-dims-netin")){
			expectReject(name, "required: dims=", input, output, null);
		}else if(name.equals("valin-vfraction")){
			expectReject(name, "Use either valin= or vfraction>0", input, output, netIn,
				"valin="+(input), "vfraction=0.25");
		}else if(name.equals("threads-simd")){
			expectReject(name, "threads>1 currently requires simd=t", input, output, netIn,
				"threads=2", "simd=f");
		}else if(name.equals("vsimd-simd")){
			expectReject(name, "vsimd=t requires simd=t", input, output, netIn,
				"vsimd=t", "simd=f");
		}else if(name.equals("epochs-zero")){
			expectReject(name, "epochs must be positive: 0", input, output, netIn, "epochs=0");
		}else if(name.equals("startepoch-negative")){
			expectReject(name, "startepoch must be >=0", input, output, netIn, "startepoch=-1");
		}else if(name.equals("printevery-negative")){
			expectReject(name, "printevery must be >=0", input, output, netIn, "printevery=-1");
		}else if(name.equals("batch-zero")){
			expectReject(name, "batch must be positive: 0", input, output, netIn, "batch=0");
		}else if(name.equals("density-zero")){
			expectReject(name, "density must be in (0,1]: 0.0", input, output, netIn, "density=0");
		}else if(name.equals("density-two")){
			expectReject(name, "density must be in (0,1]: 2.0", input, output, netIn, "density=2");
		}else if(name.equals("density1-two")){
			expectReject(name, "density1 must be <=1: 2.0", input, output, netIn, "density1=2");
		}else if(name.equals("edgeblock-zero")){
			expectReject(name, "edgeblock must be positive: 0", input, output, netIn, "edgeblock=0");
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
	}

	/** Asserts that construction rejects an invalid argument combination. */
	private static void expectReject(String name, String message, Path input, Path output,
			Path netIn, String... extra){
		try{
			construct(input, output, netIn, extra);
			throw new AssertionError(name+" should be rejected by the constructor");
		}catch(IllegalArgumentException expected){
			if(!expected.getMessage().contains(message)){throw expected;}
		}
	}

	/** Constructs a trainer with omitted dims= and the supplied extra arguments. */
	private static void construct(Path input, Path output, Path netIn, String... extra){
		final String[] args=new String[3+(netIn==null ? 0 : 1)+extra.length];
		int i=0;
		args[i++]="in="+input;
		args[i++]="out="+output;
		if(netIn!=null){args[i++]="netin="+netIn;}
		args[i++]="overwrite=t";
		for(String s : extra){args[i++]=s;}
		new RegressionTrainer(args);
	}

	/** One tiny existing input for constructor file checks; it is never parsed. */
	private static final String ROWS="#dims\t1\t1\n0\t0\n";

	private RegressionTrainerConstructorValidationTest(){}
}
