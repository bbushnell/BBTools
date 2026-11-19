package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import ddl.DDLFamilyHoldout;
import prot.DigestSuffix;
import prot.FamilyCalibration.MemberInput;
import prot.FamilyCalibration.MemberInputs;
import shared.Shared;

/** Canonical multi-bucket input, raw-consensus preservation and fail-closed preparation checks.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyPreparationInputTest{

	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Require one fresh fixture directory");}
		final Path root=Paths.get(args[0]).toAbsolutePath().normalize(); Files.createDirectory(root);
		Shared.AMINO_IN=true;
		final Fixture f=new Fixture(root);
		final MemberInputs inputs=f.open();
		check(inputs.familyCount()==3 && inputs.members==3 && inputs.buckets==2 && inputs.batches==2, "Pinned input census");
		reject(()->inputs.verifyComplete());
		final MemberInput[] left=inputs.loadBucket(0, 3, 4096);
		check(left.length==2 && left[0].rank==0 && left[1].rank==2, "Bucket preserves physical family ranks including an empty family");
		check(left[0].id.equals("f0") && left[0].members.length==2 && Arrays.equals(left[0].members[0], f.a) &&
			Arrays.equals(left[0].members[1], f.w), "Members preserve canonical order across batches");
		check(Arrays.equals(left[0].consensus, f.raw), "Original B/Z/J/X consensus bytes must not be reconstructed from collapsed encoding");
		check(left[1].members.length==0 && Arrays.equals(left[1].consensus, new byte[]{'X'}), "Empty family retains explicit input without invented members");
		reject(()->inputs.loadBucket(0, 3, 4096));
		final MemberInput[] right=inputs.loadBucket(1, 3, 4096);
		check(right.length==1 && right[0].rank==1 && right[0].members.length==1 && Arrays.equals(right[0].members[0], f.c), "Second bucket conserves its members");
		check(Arrays.equals(left[0].members[0], f.a) && Arrays.equals(left[0].members[1], f.w) && Arrays.equals(left[0].consensus, f.raw),
			"Returned arrays survive later bucket loads without scratch aliasing");
		inputs.verifyComplete();
		final MemberInputs memberCap=f.open();
		reject(()->memberCap.loadBucket(0, 1, 4096));
		reject(()->memberCap.loadBucket(1, 3, 4096));
		reject(()->memberCap.verifyComplete());
		final MemberInputs byteCap=f.open();
		reject(()->byteCap.loadBucket(0, 3, 256));
		final MemberInputs changed=f.open();
		final byte[] original=Files.readAllBytes(f.parts[0]);
		write(f.parts[0], new String(original, StandardCharsets.US_ASCII).replace('A', 'W'));
		reject(()->changed.loadBucket(0, 3, 4096));
		Files.write(f.parts[0], original);
		reject(()->changed.loadBucket(1, 3, 4096));
		final byte[] originalConsensus=Files.readAllBytes(f.consensus);
		write(f.consensus, new String(originalConsensus, StandardCharsets.US_ASCII).replace('A', 'W'));
		reject(()->inputs.verifyComplete());
		Files.write(f.consensus, originalConsensus);
		reject(()->new MemberInputs(f.families, "00000000000000000000", f.chunks, hash(f.chunks), f.consensus, hash(f.consensus), 3, 3, 2, 2, 1000));
		final String marker="FAMILY_PREPARATION_INPUT_TEST_PASS checks="+checks;
		write(root.resolve("PASS"), marker+"\n");
		System.out.println(marker);
	}

	private static final class Fixture{
		Fixture(Path root) throws Exception{
			families=root.resolve("families.tsv"); chunks=root.resolve("chunks.tsv"); consensus=root.resolve("consensus.faa");
			a=sequence('A'); w=sequence('W'); c=sequence('C'); raw=a.clone();
			System.arraycopy("BZJX".getBytes(StandardCharsets.US_ASCII), 0, raw, 96, 4);
			write(families, DDLFamilyHoldout.FAMILY_HEADER+"\n0\tf0\t2\n1\tf1\t1\n2\tempty\t0\n");
			write(consensus, ">f0\n"+ascii(raw)+"\n>f1\n"+ascii(c)+"\n>empty\nX\n");
			final String[] records={"0\tf0\tg0|a\tg0\t0\t"+ascii(a)+"\n", "1\tf1\tg0|c\tg0\t1\t"+ascii(c)+"\n",
				"0\tf0\tg1|w\tg1\t0\t"+ascii(w)+"\n", ""};
			final StringBuilder manifest=new StringBuilder(DDLFamilyHoldout.CHUNK_HEADER).append('\n');
			for(int i=0; i<4; i++){
				parts[i]=root.resolve("part"+i+".tsv"); write(parts[i], DDLFamilyHoldout.MEMBER_HEADER+"\n"+records[i]);
				final int count=i==3 ? 0 : 1;
				manifest.append(String.format(Locale.ROOT, "b%02d", i/2)).append('\t').append(i%2).append('\t').append(parts[i])
					.append('\t').append(hash(parts[i])).append('\t').append(count).append('\t').append(count*100).append('\n');
			}
			write(chunks, manifest.toString());
		}
		MemberInputs open(){return new MemberInputs(families, hash(families), chunks, hash(chunks), consensus, hash(consensus), 3, 3, 2, 2, 1000);}
		final Path families, chunks, consensus;
		final Path[] parts=new Path[4];
		final byte[] a, w, c, raw;
	}
	private static byte[] sequence(char symbol){final byte[] result=new byte[100]; Arrays.fill(result, (byte)symbol); return result;}
	private static String ascii(byte[] bytes){return new String(bytes, StandardCharsets.US_ASCII);}
	private static String hash(Path path){return DigestSuffix.file(path.toString());}
	private static void write(Path path, String text) throws Exception{Files.write(path, text.getBytes(StandardCharsets.US_ASCII));}
	private static void reject(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){checks++; return;}
		throw new AssertionError("Invalid preparation input/session was accepted");
	}
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);} checks++;}
	private static int checks;
}
