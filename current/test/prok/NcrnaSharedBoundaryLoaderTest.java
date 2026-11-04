package test.prok;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Locale;

import dna.Data;
import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.NcrnaBoundaryResources;
import prok.NcrnaBoundaryScorer;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;
import prok.TrnaBoundaryFeatures.BoundaryType;
import prok.TrnaNinemerTableBuilder;
import prok.NcrnaTextIO.Input;
import shared.Shared;

/** Offline shared-boundary checks using current packaged resources and literal lengths.
 * No retired installation, downloaded input, trained fixture or resource mutation is needed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSharedBoundaryLoaderTest{

	/*--------------------------------------------------------------*/
	/*----------------          Entry Point         ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSharedBoundaryLoaderTest.class, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			final Parser parser=new Parser();
			String selected=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("family")){selected=b;}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown loader-fixture option: "+arg);}
			}
			Shared.setThreads(1);
			final Path scratch=Files.createTempDirectory("ncrna-shared-loader-");
			try(AutoCloseable cleanup=()->Files.delete(scratch)){
				int tested=0;
				for(final String id:FAMILIES){
					if(selected==null || selected.equals(id)){
						checkFamily(id, scratch);
						tested++;
					}
				}
				check(tested>0, "Select rnasep, tmrna, srp_small or srp_large");
				System.out.println("Shared boundary loader passed: families="+tested+
					"; current-package binding, literal median, table types/counts, predictions, overrides and disabled loading");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Checks            ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercises the public loader and adapter against their packaged payloads. */
	private static void checkFamily(final String id, final Path scratch){
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final NcrnaFamilyLoader.Options ordinary=new NcrnaFamilyLoader.Options();
		final String manifest=Data.findPath("?"+NcrnaFamilyRegistry.get(id).directory+"/manifest.tsv", false);
		check(manifest!=null, "Current packaged manifest is required for the loader fixture: "+id);
		ordinary.resources.root=new File(manifest).getParent();
		ordinary.boundaryOverride=false;
		final NcrnaFamily family=NcrnaFamilyLoader.load(id, ordinary, resolver);
		check(family.library.length>0 && family.models.length==family.library.length,
			"Installed consensus and HBM axes must remain nonempty and aligned: "+id);
		check(family.boundary5NetTemplate==null && family.boundary3NetTemplate==null,
			"Explicit boundary=false must not load boundary networks: "+id);
		final String network=resolver.resolve(id, "boundary", ordinary.resources).path;
		final String start=resolver.resolve(id, "boundary_start", ordinary.resources).path;
		final String stop=resolver.resolve(id, "boundary_stop", ordinary.resources).path;
		final NcrnaBoundaryResources loaded=new NcrnaBoundaryResources(id, family.library, ordinary.resources, resolver);
		final CellNet direct=NcrnaBoundaryScorer.load(network);
		check(direct.numInputs()==10 && loaded.net.numInputs()==10 && direct.numOutputs()==1
			&& loaded.net.numOutputs()==1, "Shared boundary inference requires ten inputs and one output: "+id);
		final float[] input=new float[10];
		for(int row=0; row<8; row++){
			for(int i=0; i<input.length; i++){input[i]=(float)Math.sin((i+1)*(row+1));}
			direct.applyInput(input);
			loaded.net.applyInput(input);
			direct.feedForward();
			loaded.net.feedForward();
			check(Arrays.equals(direct.getOutput(), loaded.net.getOutput()) && Float.isFinite(loaded.net.getOutput()[0]),
				"The adapter must preserve the current network's finite predictions: "+id);
		}
		checkTable(start, loaded.start, BoundaryType.START);
		checkTable(stop, loaded.stop, BoundaryType.STOP);
		final byte[][] lengths={new byte[10], new byte[20], new byte[30], new byte[100]};
		check(new NcrnaBoundaryResources(id, lengths, ordinary.resources, resolver).meanLen==30,
			"The upper median of literal lengths 10, 20, 30, 100 is 30, not 20 or 40");
		final FamilyResourceResolver.Overrides swapped=new FamilyResourceResolver.Overrides();
		swapped.root=ordinary.resources.root;
		swapped.boundaryStart=stop;
		fails(()->new NcrnaBoundaryResources(id, family.library, swapped, resolver), "start table");
		final FamilyResourceResolver.Overrides missing=new FamilyResourceResolver.Overrides();
		missing.root=ordinary.resources.root;
		missing.boundaryStart=scratch.resolve("absent.tsv").toString();
		fails(()->new NcrnaBoundaryResources(id, family.library, missing, resolver), "boundary_start");

		final NcrnaFamilyLoader.Options enabled=new NcrnaFamilyLoader.Options();
		enabled.resources.root=ordinary.resources.root;
		enabled.boundary=true;
		enabled.boundaryOverride=true;
		final NcrnaFamily refined=NcrnaFamilyLoader.load(id, enabled, resolver);
		check(refined.boundary5NetTemplate!=null && refined.boundary3NetTemplate!=null,
			"Explicit boundary loading must bind both ends: "+id);
		check(Arrays.deepEquals(family.library, refined.library) && Arrays.equals(family.modelNames, refined.modelNames),
			"Boundary selection must not replace or reorder the core model library: "+id);
		final NcrnaFamilyLoader.Options disabled=new NcrnaFamilyLoader.Options();
		disabled.enabled=false;
		disabled.resources.root=scratch.resolve("absent").toString();
		check(NcrnaFamilyLoader.load(id, disabled, resolver)==null, "Disabled families must not touch an absent resource root");
		disabled.enabled=true;
		fails(()->NcrnaFamilyLoader.load(id, disabled, resolver), "manifest");
		final String[] names=family.modelNames.clone();
		resolver.resolve(id.equals("rnasep") ? "tmrna" : "rnasep", "seeds", null);
		check(start.equals(resolver.resolve(id, "boundary_start", ordinary.resources).path) && Arrays.equals(names, family.modelNames),
			"Resolver family switching must restore the requested binding without mutating an existing library: "+id);
	}

	/** Checks every stored count against the loaded table, with the declared end checked independently. */
	private static void checkTable(final String path, final TrnaNinemerTableBuilder.LoadedTable loaded, final BoundaryType type){
		check(loaded.type==type && loaded.k>0, "Boundary table must match its declared endpoint: "+path);
		final LineParser1 parser=new LineParser1('\t');
		int rows=0;
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				parser.set(line);
				if(parser.terms()!=2 || parser.setBounds(0)!=loaded.k){
					throw new AssertionError("Expected one complete k-mer and count at row "+rows+" in "+path);
				}
				final int count=parser.parseInt(1);
				if(loaded.table.rawCount(line, 0)!=count){
					throw new AssertionError("Loaded count differs from stored row "+rows+" in "+path);
				}
				rows++;
			}
		}
		check(rows>0, "The current packaged table must contain actual count rows: "+path);
	}

	private static void fails(final Runnable action, final String message){
		try{action.run();}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(message), "Wrong rejection: "+expected);
			return;
		}
		throw new AssertionError("Expected required-resource rejection: "+message);
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] FAMILIES={"rnasep", "tmrna", "srp_small", "srp_large"};
}
