package org.panda.resource;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A tool to get information on drugs and their targets.
 * @author Ozgun Babur
 */
public class PiHelper extends FileServer
{
	private Map<String, Set<String>> id2Target;
	private Map<String, Set<String>> target2id;
	private Map<String, String> id2Name;
	private Map<String, String> name2id;

	private Map<String, String> description;
	private Set<String> fdaApproved;
	private Set<String> cancerDrug;
	private Set<String> nutraceutical;

	private static PiHelper instance;
	public static PiHelper get()
	{
		if (instance == null) instance = new PiHelper();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"drugs.tsv", "drugtargets.tsv"};
	}

	@Override
	public boolean load() throws IOException
	{
		id2Name = new HashMap<>();
		name2id = new HashMap<>();
		description = new HashMap<>();
		fdaApproved = new HashSet<>();
		cancerDrug = new HashSet<>();
		nutraceutical = new HashSet<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(4).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			String name = t[1];

			id2Name.put(id, name);
			name2id.put(name, id);
			name2id.put(name.toLowerCase(), id);

			for (String syn : t[2].split(";"))
			{
				name2id.put(syn, id);
				name2id.put(syn.toLowerCase(), id);

				int a = syn.indexOf(" (");
				if (a > 0)
				{
					int b = syn.indexOf(")");

					if (b > a)
					{
						String s = syn.substring(0, a).trim();
						name2id.put(s, id);
						name2id.put(s.toLowerCase(), id);
					}
				}
			}

			description.put(id, t[3]);

			if (Boolean.valueOf(t[6])) fdaApproved.add(id);
			if (Boolean.valueOf(t[7])) cancerDrug.add(id);
			if (Boolean.valueOf(t[8])) nutraceutical.add(id);
		});

		id2Target = new HashMap<>();
		target2id = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[1]).skip(4).map(l -> l.split("\t")).forEach(t ->
		{
			String id = name2id.get(t[2]);
			String target = t[1];
			if (!id2Target.containsKey(id)) id2Target.put(id, new HashSet<>());
			id2Target.get(id).add(target);
			if (!target2id.containsKey(target)) target2id.put(target, new HashSet<>());
			target2id.get(target).add(id);
		});

		return true;
	}

	void printGenesTargetedByCancerDrugs() throws IOException
	{
		Map<String, Set<String>> gene2drug = new HashMap<>();
		Set<String> set = Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/List-all-drugs.txt")).skip(2)
			.map(l -> name2id.get(l.split("\t")[1].trim())).filter(Objects::nonNull).collect(Collectors.toSet());
		set.addAll(cancerDrug);
		for (String id : set)
		{
			String name = id2Name.get(id);
			if (id2Target.containsKey(id))
			{
				for (String gene : id2Target.get(id))
				{
					if (!gene2drug.containsKey(gene)) gene2drug.put(gene, new HashSet<>());
					gene2drug.get(gene).add(name);
				}
			}
		}
		gene2drug.keySet().stream().sorted().forEach(gene -> System.out.println(gene + "\t" + gene2drug.get(gene)));
	}

	void checkExistenceOfDrugs(Set<String> set)
	{
		for (String name : set)
		{
			if (!name2id.containsKey(name.toLowerCase()))
			{
				System.out.println("Drug not found = " + name);
			}
		}
	}

	public static void main(String[] args) throws IOException
	{
		get().printGenesTargetedByCancerDrugs();
//		Set<String> set = Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/List-all-drugs.txt"))
//			.skip(2).map(l -> l.split("\t")[1].trim()).collect(Collectors.toSet());
//		get().checkExistenceOfDrugs(set);
	}
}
