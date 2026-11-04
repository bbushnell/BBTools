package test.prok;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.NcrnaFamily;
import prok.NcrnaFamilyExport;
import prok.NcrnaFamilyExport.Boundary;
import prok.NcrnaFamilyLoader;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Finite export controls using installed shared-boundary libraries.
 * Tests reread output through the actual consumer; the production exporter does not.
 * Output remains in the explicit fresh fixture directory for independent review.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaFamilyExportTest {

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, null, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}

	/** Main owns the selected log across every fixture path, including early returns. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		boolean nestedOnly=false;
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(java.util.Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("nestedonly")){nestedOnly=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		if(parser.out1==null){throw new IllegalArgumentException("Specify a fresh out= fixture directory");}
		Shared.setThreads(1);
		ReadWrite.setZipThreads(1);
		final Path root=Files.createDirectory(new File(parser.out1).toPath());
		checkNestedConfiguration(root);
		if(nestedOnly){
			System.out.println("Nested sibling configuration passed: config-relative, parent-relative, wrong binding, missing binding and absolute dependency");
			return;
		}
		checkOff(root, "fmn", "", null);
		checkOff(root, "sixs_rf01685", "RF01685", "sixs_rf00013");
		checkOff(root, "srp_small", "small", "srp_large");
		checkOff(root, "cobalamin_rf00174", "RF00174", "cobalamin_rf01482");
		checkOn(root);
		checkChangedSharedSeeds(root);
		System.out.println("Family exports passed: ungrouped/grouped OFF, SAM-I ON, native loader geometry, tiny fixed decimals, gzip, sibling payloads and shared-seed rejection");
	}

	private static void checkOff(final Path root, final String id, final String member, final String sibling) throws IOException{
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final String consensus=resolver.resolve(id, "consensus", null).path;
		final String models=resolver.resolve(id, "models", null).path;
		final String seeds=resolver.resolve(id, "seeds", null).path;
		// Different compressed representation, same forward seed membership.
		final Path compressed=root.resolve(id+".seeds.fa.gz");
		copyText(seeds, compressed.toString());
		final NcrnaFamilyLoader.Options originalOptions=options(null);
		originalOptions.resources.seeds=compressed.toString();
		final NcrnaFamily original=NcrnaFamilyLoader.load(id, originalOptions, resolver);
		final Path destination=root.resolve(id+" exported");
		final Boundary policy=new Boundary(null, null, null, 2, 3, Float.MIN_VALUE, .15f, Float.NaN);
		NcrnaFamilyExport.write(original, consensus, models, compressed.toString(), policy, destination.toFile());
		final NcrnaFamilyLoader.Options loadedOptions=options(destination);
		loadedOptions.boundary=true;
		final NcrnaFamily loaded=NcrnaFamilyLoader.load(id, loadedOptions, new FamilyResourceResolver());
		assert(loaded.boundary5NetTemplate==null && loaded.boundary3NetTemplate==null)
			: "Exported OFF must replace the loader's requested default without loading networks";
		assert(loaded.boundaryStartOffsets.length==5 && loaded.boundaryStopOffsets.length==7
			&& loaded.boundaryMarginStart==Float.MIN_VALUE && loaded.boundaryMarginStop==.15f)
			: "Native reload must preserve declared geometry and the smallest positive float margin";
		assert(Arrays.deepEquals(original.library, loaded.library) && Arrays.equals(original.modelNames, loaded.modelNames)
			&& original.kmerSet.size()==loaded.kmerSet.size()) : "Core library and complete seed population must survive export";
		for(long key:original.kmerSet.toArray()){assert(loaded.kmerSet.contains(key)) : "Every accepted seed must survive export";}
		final Path selected=member.isEmpty() ? destination : destination.resolve(member);
		assert(!Files.exists(selected.resolve("boundary.bbnet")) && !Files.exists(selected.resolve("tables")))
			: "OFF exports must not accidentally retain old selected-family boundary products";
		checkMetadata(destination.resolve("manifest.tsv"));
		if(sibling!=null){
			final NcrnaFamily before=NcrnaFamilyLoader.load(sibling, options(null), new FamilyResourceResolver());
			final NcrnaFamily after=NcrnaFamilyLoader.load(sibling, options(destination), new FamilyResourceResolver());
			assert(Arrays.deepEquals(before.library, after.library) && Arrays.equals(before.modelNames, after.modelNames))
				: "Unselected sibling models must remain unchanged";
			final Path installed=new File(consensus).getAbsoluteFile().getParentFile().getParentFile().toPath();
			checkPreservedBindings(installed, destination, member);
		}
		boolean rejected=false;
		try{NcrnaFamilyExport.write(original, consensus, models, seeds, policy, destination.toFile());}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("fresh output");}
		assert(rejected) : "Existing output must never be overwritten by an export";
	}

	private static void checkOn(final Path root) throws IOException{
		final String id="sam1";
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final NcrnaFamily original=NcrnaFamilyLoader.load(id, options(null), resolver);
		final Path net=root.resolve("test network.bbnet.gz");
		copyText(resolver.resolve(id, "boundary", null).path, net.toString());
		final String start=resolver.resolve(id, "boundary_start", null).path;
		final String stop=resolver.resolve(id, "boundary_stop", null).path;
		final Path destination=root.resolve("sam1 on");
		final Boundary policy=new Boundary(net.toString(), start, stop, 8, 12, .125f, .25f, Float.MIN_VALUE);
		NcrnaFamilyExport.write(original, resolver.resolve(id, "consensus", null).path,
			resolver.resolve(id, "models", null).path, resolver.resolve(id, "seeds", null).path, policy, destination.toFile());
		final NcrnaFamily loaded=NcrnaFamilyLoader.load(id, options(destination), new FamilyResourceResolver());
		assert(loaded.boundary5NetTemplate!=null && loaded.boundary3NetTemplate!=null
			&& loaded.boundaryStartOffsets.length==17 && loaded.boundaryStopOffsets.length==25
			&& loaded.boundaryMarginStart==.125f && loaded.boundaryMarginStop==.25f)
			: "Native ON reload must preserve the gzip network, tables, radii and margins";
		assert(Files.isRegularFile(destination.resolve("boundary.bbnet.gz"))) : "Copied gzip network must retain compression dispatch";
		checkMetadata(destination.resolve("manifest.tsv"));
		try(Input input=new Input(destination.resolve("boundary.config").toString())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				final String value=new String(line, java.nio.charset.StandardCharsets.UTF_8);
				if(value.startsWith("cutoff=")){
					assert(!value.contains("E") && !value.contains("e") && Float.parseFloat(value.substring(7))==Float.MIN_VALUE)
						: "Numeric output must be fixed decimal and preserve the subnormal cutoff";
				}
			}
		}
	}

	private static void checkChangedSharedSeeds(final Path root){
		final String id="srp_small";
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final Path seeds=root.resolve("different.fa");
		try(Output output=new Output(seeds.toString())){output.writer.print(new ByteBuilder(">different\nAAAAAAAAAAAAAAAAA\n"));}
		final NcrnaFamilyLoader.Options options=options(null);
		options.resources.seeds=seeds.toString();
		final NcrnaFamily original=NcrnaFamilyLoader.load(id, options, resolver);
		final Path destination=root.resolve("must not exist");
		boolean rejected=false;
		try{
			NcrnaFamilyExport.write(original, resolver.resolve(id, "consensus", null).path,
				resolver.resolve(id, "models", null).path, seeds.toString(),
				new Boundary(null, null, null, 0, 0, 0f, 0f, Float.NaN), destination.toFile());
		}catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("shared seeds");}
		assert(rejected && !Files.exists(destination)) : "Changed shared seeds must fail before creating a partial export";
	}

	/** No installed resource mutation: literal paths distinguish config-parent from member-root resolution. */
	private static void checkNestedConfiguration(final Path root) throws IOException{
		final Path member=Files.createDirectory(root.resolve("nested member"));
		final Path configs=Files.createDirectory(member.resolve("configs"));
		final Path config=configs.resolve("boundary.config");
		final Path nestedNet=configs.resolve("boundary.bbnet");
		final Path rootNet=member.resolve("boundary.bbnet");
		Files.write(nestedNet, new byte[]{1});
		Files.write(rootNet, new byte[]{2});
		writeConfig(config, "  # valid whitespace comment\nnet=boundary.bbnet\n");
		NcrnaFamilyExport.requirePortableConfiguration(config.toString(), new String[]{null, nestedNet.toString()});
		rejectConfig(config, new String[]{rootNet.toString()}, "not manifest-bound");
		rejectConfig(config, new String[0], "not manifest-bound");
		Files.delete(config);
		writeConfig(config, "net=../boundary.bbnet\n");
		NcrnaFamilyExport.requirePortableConfiguration(config.toString(), new String[]{rootNet.toString()});
		rejectConfig(config, new String[]{nestedNet.toString()}, "not manifest-bound");
		Files.delete(config);
		writeConfig(config, "net="+nestedNet.toAbsolutePath()+"\n");
		rejectConfig(config, new String[]{nestedNet.toString()}, "must be relative");
	}

	private static void writeConfig(final Path path, final String text){
		try(Output output=new Output(path.toString())){output.writer.print(new ByteBuilder(text));}
	}

	private static void rejectConfig(final Path path, final String[] bindings, final String message){
		boolean rejected=false;
		try{NcrnaFamilyExport.requirePortableConfiguration(path.toString(), bindings);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(message);}
		assert(rejected) : "Wrong/external sibling dependencies must be rejected by the intended path-binding condition";
	}

	/** Independent raw-byte comparison of every unselected/shared manifest payload. */
	private static void checkPreservedBindings(final Path original, final Path exported, final String selected) throws IOException{
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(original.resolve("manifest.tsv").toString())){
			input.reader.nextLine();
			input.reader.nextLine();
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				assert(parser.terms()==8) : "Group fixture must use the actual eight-column manifest";
				final String member=parser.parseString(0);
				if(member.equals(selected)){continue;}
				final String relative=parser.parseString(4);
				final Path a=(member.equals(".") ? original : original.resolve(member)).resolve(relative);
				final Path b=(member.equals(".") ? exported : exported.resolve(member)).resolve(relative);
				assert(sameBytes(a, b)) : "An unselected/shared payload changed: "+member+"/"+relative;
			}
		}
	}

	private static boolean sameBytes(Path a, Path b) throws IOException{
		if(Files.size(a)!=Files.size(b)){return false;}
		try(InputStream first=Files.newInputStream(a); InputStream second=Files.newInputStream(b)){
			final byte[] x=new byte[65536];
			final byte[] y=new byte[65536];
			for(int count=first.read(x); count>=0; count=first.read(x)){
				int offset=0;
				while(offset<count){
					final int n=second.read(y, offset, count-offset);
					if(n<0){return false;}
					offset+=n;
				}
				for(int i=0; i<count; i++){if(x[i]!=y[i]){return false;}}
			}
			return second.read()<0;
		}
	}

	private static void checkMetadata(final Path path){
		try(Input input=new Input(path.toString())){
			input.reader.nextLine();
			input.reader.nextLine();
			final LineParser1 parser=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				assert(parser.termEquals("-", parser.terms()-2) && parser.termEquals("-", parser.terms()-1))
					: "New manifest metadata must remain absent rather than hashed or counted";
			}
		}
	}

	private static void copyText(String source, String target){
		try(Input input=new Input(source); Output output=new Output(target)){
			final ByteBuilder text=new ByteBuilder();
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				text.clear();
				output.writer.print(text.append(line).nl());
			}
		}
	}

	private static NcrnaFamilyLoader.Options options(Path root){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.endpoints=false;
		if(root!=null){options.resources.root=root.toString();}
		return options;
	}

	private static final class Input implements AutoCloseable {
		Input(String path){reader=ByteFile.makeByteFile(path, false);}
		@Override
		public void close(){if(reader.close()){throw new IllegalStateException("Fixture input close failed");}}
		final ByteFile reader;
	}
}
