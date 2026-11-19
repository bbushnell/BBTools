package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Locale;

import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CallGenes;
import prok.ConservedRnaSeedIndex;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyExport;
import prok.NcrnaFamilyLoader;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Runs the ordinary loader with a different seed set for exactly one B12 member.
 * Each invocation uses a fresh JVM. Both the explicit seed and explicit member-root
 * routes must preserve the two untouched shipped sets and actual shared-index hits.
 * No genome benchmark, model training or scientific seed selection is performed.
 * Adapted from template.A_SampleBasic2 and CaptureNcrnaWindowsTest.
 * @author Brian Bushnell, Raiden
 */
public final class CallGenesCobalaminSeedOverrideTest{

	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, CallGenesCobalaminSeedOverrideTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String member=null;
			boolean rootOverride=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("member")){member=CallGenes.parseNcrnaFamily(b);}
				else if(a.equals("root")){rootOverride=Parse.parseBoolean(b);}
				else{check(parser.parse(arg, a, b), "Unknown scoped-seed fixture argument: "+arg);}
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory member=cobalamin_rfNNNNN root=t|f");
			new CallGenesCobalaminSeedOverrideTest(new File(parser.out1).getAbsoluteFile(), member, rootOverride).run();
			pp.outstream.println("CALLGENES_COBALAMIN_SEED_OVERRIDE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private CallGenesCobalaminSeedOverrideTest(final File root_, final String member_, final boolean rootOverride_){
		root=root_;
		member=member_;
		rootOverride=rootOverride_;
		for(int i=0; i<NAMES.length; i++){if(NAMES[i].equals(member)){selected=i;}}
		check(selected>=0 && root.mkdir(), "A fresh fixture and exactly one B12 member are required");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Ordinary Loader Test    ----------------*/
	/*--------------------------------------------------------------*/

	private void run() throws IOException{
		check(GeneCaller.ncrnaFamilies.isEmpty(), "Run each scoped-seed case in a fresh JVM before loading caller resources");
		final NcrnaFamily[] shipped=new NcrnaFamily[NAMES.length];
		for(int i=0; i<NAMES.length; i++){shipped[i]=NcrnaFamilyLoader.load(NAMES[i]);}
		CallGenes.requireCobalaminGroup(shipped);
		final long[] replacement=new long[2];
		int found=0;
		for(long key=0; found<replacement.length && key<1024; key++){
			if(!shipped[0].kmerSet.contains(key) && !shipped[1].kmerSet.contains(key) && !shipped[2].kmerSet.contains(key)){
				replacement[found++]=key;
			}
		}
		check(found==replacement.length, "The fixture needs two keys absent from every original seed set");
		final String override;
		if(rootOverride){
			final FamilyResourceResolver resolver=new FamilyResourceResolver();
			final String consensus=resolver.resolve(member, "consensus", null).path;
			final String models=resolver.resolve(member, "models", null).path;
			final String seeds=resolver.resolve(member, "seeds", null).path;
			final File exported=new File(root, "member-root");
			NcrnaFamilyExport.write(shipped[selected], consensus, models, seeds,
				new NcrnaFamilyExport.Boundary(null, null, null, 2, 2, 0, 0, Float.NaN), exported);
			final File shared=new File(exported, new File(seeds).getName());
			check(shared.isFile(), "The grouped export must retain the shared seed payload before the scoped replacement");
			Files.move(shared.toPath(), new File(exported, "original-seeds.saved").toPath());
			writeSeeds(shared.getPath(), replacement);
			override="ncrnafamilyroot="+exported.getPath();
		}else{
			writeSeeds(path("replacement.fa"), replacement);
			override="ncrnakmers="+path("replacement.fa");
		}
		write(path("input.fa"), new ByteBuilder(">no_seed\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"));
		CallGenes.main(new String[]{"in="+path("input.fa"), "out=null", "cdsonly=t", "callcds=f",
			"ncrna=f", "cobalamin=t", "tmrna=f", "sixs=f", "fmn=f", "passes=1", "taxonomy=f",
			"ncrnafamily="+member, override, "t=1"});
		final NcrnaFamily[] loaded=new NcrnaFamily[NAMES.length];
		int members=0;
		for(final NcrnaFamily family:GeneCaller.ncrnaFamilies){
			for(int i=0; i<NAMES.length; i++){
				if(!NAMES[i].equals(family.name)){continue;}
				check(loaded[i]==null, "The ordinary caller must not publish a B12 member twice: "+family.name);
				loaded[i]=family;
				members++;
			}
		}
		check(members==3, "The ordinary combined caller must retain all three named B12 members alongside other ncRNAs");
		for(int i=0; i<loaded.length; i++){
			check(loaded[i].name.equals(NAMES[i]) && loaded[i].kLong==17 && loaded[i].library.length==shipped[i].library.length,
				"Scoped seed replacement must preserve every identity, seed length and model roster");
			final LongHashSet expected=i==selected ? replacementSet(replacement) : shipped[i].kmerSet;
			equalSeeds(expected, loaded[i].kmerSet);
		}
		final ConservedRnaSeedIndex index=new ConservedRnaSeedIndex(
			new LongHashSet[]{loaded[0].kmerSet, loaded[1].kmerSet, loaded[2].kmerSet});
		final ByteBuilder sequence=new ByteBuilder();
		for(final long key:replacement){
			sequence.clear().appendKmer(key, 17);
			final ConservedRnaSeedIndex.ScanResult scan=index.scan(sequence.toBytes());
			for(int i=0; i<loaded.length; i++){
				check(scan.hits(i).length==(i==selected ? 1 : 0),
					"A replacement key must reach exactly its selected family slot in the actual shared index");
			}
		}
		sequence.clear().appendKmer(shipped[0].kmerSet.toArray()[0], 17);
		final ConservedRnaSeedIndex.ScanResult retained=index.scan(sequence.toBytes());
		for(int i=0; i<loaded.length; i++){
			check(retained.hits(i).length==(i==selected ? 0 : 1),
				"An original shared key must still reach both untouched siblings after a scoped replacement");
		}
		boolean rejected=false;
		try{CallGenes.requireCobalaminGroup(loaded);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "The default/group validator must remain strict outside the explicit single-member loading route");
	}

	private static LongHashSet replacementSet(final long[] keys){
		final LongHashSet set=new LongHashSet(4);
		for(final long key:keys){set.add(key);}
		return set;
	}
	private static void equalSeeds(final LongHashSet expected, final LongHashSet actual){
		check(expected.size()==actual.size(), "Scoped replacement changed an untouched seed denominator");
		for(final long key:expected.toArray()){check(actual.contains(key), "Scoped replacement silently changed an untouched seed key");}
	}
	private static void writeSeeds(final String path, final long[] keys){
		final ByteBuilder text=new ByteBuilder();
		for(int i=0; i<keys.length; i++){text.append(">replacement_").append(i).nl().appendKmer(keys[i], 17).nl();}
		write(path, text);
	}
	private static void write(final String path, final ByteBuilder text){try(Output out=new Output(path)){out.writer.print(text);}}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private final String member;
	private final boolean rootOverride;
	private int selected=-1;
	private static final String[] NAMES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
}
