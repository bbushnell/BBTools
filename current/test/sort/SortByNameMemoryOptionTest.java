package test.sort;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import sort.SortByName;

/** Constructor-only coverage for SortByName memory option validation.
 * @author Jean
 */
public final class SortByNameMemoryOptionTest{

	/** Runs baseline|fixed in a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		if(args.length!=2 || (!args[0].equals("baseline") && !args[0].equals("fixed"))){
			throw new IllegalArgumentException("Expected baseline|fixed and a fresh fixture directory");
		}
		final boolean fixed=args[0].equals("fixed");
		final Path root=Paths.get(args[1]).toAbsolutePath();
		Files.createDirectory(root);
		System.out.println("case\taccepted\texpected");
		for(Case c : CASES){
			check(root, c, fixed);
		}
	}

	/** Constructs SortByName and verifies whether this option is accepted. */
	private static void check(Path root, Case c, boolean fixed) throws Exception{
		final Path dir=root.resolve(c.name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		Files.write(input, READ.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final boolean expected=(!fixed || c.valid);
		boolean accepted=false;
		try{
			new SortByName(new String[]{"in="+input, "out="+output, "delete=f",
				"overwrite=f", c.option});
			accepted=true;
		}catch(IllegalArgumentException expectedFailure){
			if(expected){
				throw new AssertionError("Unexpected rejection for "+c.name+": "+c.option,
					expectedFailure);
			}
			final String message=expectedFailure.getMessage();
			assert(message!=null && message.startsWith(c.prefix)) :
				c.name+" rejected for the wrong reason: "+message;
		}
		assert(accepted==expected) : c.name+" accepted="+accepted+" expected="+expected;
		System.out.println(c.name+'\t'+accepted+'\t'+expected);
	}

	/** One constructor option case. */
	private static final class Case{
		Case(String name_, String option_, boolean valid_, String prefix_){
			name=name_;
			option=option_;
			valid=valid_;
			prefix=prefix_;
		}

		final String name, option;
		final boolean valid;
		final String prefix;
	}

	private static final Case[] CASES={
		new Case("negative-memlimit", "memlimit=-0.25", false, "memlimit must be finite and nonnegative"),
		new Case("nan-memlimit", "memlimit=NaN", false, "memlimit must be finite and nonnegative"),
		new Case("infinite-memlimit", "memlimit=Infinity", false, "memlimit must be finite and nonnegative"),
		new Case("zero-memlimit", "memlimit=0", true, ""),
		new Case("negative-memmult", "memmult=-0.25", false, "memmult must be a finite fraction in (0, 1]"),
		new Case("zero-memmult", "memmult=0", false, "memmult must be a finite fraction in (0, 1]"),
		new Case("nan-memmult", "memmult=NaN", false, "memmult must be a finite fraction in (0, 1]"),
		new Case("infinite-memmult", "memmult=Infinity", false, "memmult must be a finite fraction in (0, 1]"),
		new Case("above-one-memmult", "memmult=1.5", false, "memmult must be a finite fraction in (0, 1]"),
		new Case("one-memmult", "memmult=1", true, ""),
	};

	private static final String READ="@r\nACGT\n+\n!!!!\n";

}
