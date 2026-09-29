# Trabalhos relacionados — levantamento para o artigo

*Esquemas de subdivisão estelar, Maubach em malhas iniciais arbitrárias, códigos LPT, e triangulação de superfícies de Riemann em CP².*
Levantamento feito em 2026-09-28. Chaves BibTeX em `related_work.bib`. Itens marcados **[verificar]** não foram lidos no original.

---

## 0. Resumo executivo (leia primeiro)

1. **O concorrente direto é Diening–Gehring–Storn (DGS, FoCM 2025).** Eles resolvem "Maubach em qualquer malha inicial, em qualquer dimensão" com uma ideia muito próxima do nosso mecanismo de "rank global": uma **coloração própria dos vértices** com N+1 cores (N ≥ n). Com N = n, a coloração é exatamente o nosso "flag/rank global" dos Teoremas 9–10 da tese. Com N > n (sempre possível, de forma gulosa e em tempo linear), eles tratam cada n-simplexo como face de um N-simplexo virtual. O artigo precisa se posicionar com clareza em relação a DGS.

2. **Resultado novo desta pesquisa: a triangulação de Gaifullin é *balanceada*.** Os 15 vértices se dividem em 5 classes de cor com 3 vértices cada. Os vértices de uma mesma classe não são adjacentes: o complemento do 1-esqueleto é formado por 5 triângulos disjuntos. A ordem de Maubach que nossa busca encontrou coincide com essa coloração: a posição de cada vértice em todas as 108 células é uma função global do vértice, e todo par de facetas compartilhadas tem i = j (verificado: 270/270 facetas). Dois fatos confirmam isso de forma independente:
   - todo triângulo está em 4 ou 6 células (grau par). Esse é o critério de Traxler/Joswig para variedades simplesmente conexas, e CP² é simplesmente conexo;
   - o grafo dual é bipartido, que é a "chess colouring" do título do artigo de Gaifullin.

   **Consequência:** a compatibilidade de Gaifullin **não** é um fenômeno fora do mecanismo de rank global. É exatamente esse mecanismo, com o rank dado por uma coloração que existe mas não vem "de construção". O relatório `report/maubach_cp2` diz o contrário na introdução e na §3 e precisa de correção. Em compensação, isso dá ao artigo um enunciado limpo: **balanceada ⇒ compatível com Maubach, com i = j** (Traxler 1997; DGS §3.1). A busca por coloração da tese passa a ser um algoritmo para *encontrar* uma (n+1)-coloração.

3. **O caso de Kühnel também se explica:** CP²₉ é 3-vizinhal, então seu 1-esqueleto é o K₉ e precisa de 9 cores. Ele não é 5-colorível, logo não é balanceado. É o caso N = 8 > n = 4 de DGS.

4. **Os Teoremas 9 e 10 da tese não são novos como fatos.** Produto de simplexos (Kuhn/Freudenthal) e subdivisão baricêntrica são os exemplos clássicos de complexos balanceados, e "balanceado ⇒ refletido" já está em Traxler (1997). A pré-subdivisão de Kossaczký (1994, n=3) e Stevenson (2008, App. A) é essencialmente a subdivisão baricêntrica sem o nível dos pontos médios de aresta: (n+1)!/2 simplexos contra (n+1)!. O que a tese tem de próprio é **o formalismo**: esquema = invariante combinatório + algoritmo, com prova indutiva nível a nível (Lemas 5–7, Teorema 8), unificando Mitchell (faceta) e Maubach (aresta) numa só linguagem estelar. É esse o ângulo a vender.

5. **O lado LPT tem análogos fortes** (t8code/TM-index, Concurrent Binary Trees, diamantes), mas nenhum combina: bisseção **conforme**, dimensão arbitrária, **semente não cúbica**, sem ponteiros. O argumento limpo para o artigo é: *se a semente é balanceada, cada célula-semente é uma cópia afim do simplexo de Kuhn de referência, com rótulos coerentes, e a travessia entre sementes preserva o índice local da faceta.* Por isso o glpt só precisa de `seed_index` e de uma tabela de adjacência da semente.

6. **Na aplicação**, a lacuna real que ocupamos é triangular a curva **intrinsecamente em CP²** (compacto, sem cartas, sem cortes de ramificação, sem tratamento especial do infinito), por continuação combinatória numa malha de CP² refinada por Maubach. As comparações naturais são:
   - Weigle–Banks (1996): contornos de funções complexas em grade de R⁴;
   - Allgower–Schmidt e Boissonnat–Kachanovich–Wintraecken: PL-continuação em Freudenthal–Kuhn;
   - Nieser–Poelke–Polthier (2010): malhas de Riemann a partir de pontos de ramificação;
   - Bertini_real: decomposição celular numérica.

   Nosso ponto fraco declarado: **não há certificação topológica** (taxa de células ruins > 0). Plantinga–Vegter, Mourrain–Pavone e Boissonnat–Wintraecken são as referências para dizer isso com honestidade.

---

## 1. Refinamento por bisseção: newest vertex / Maubach

### 1.1 Linhagem clássica

| Trabalho | O que faz | Condição na malha inicial |
|---|---|---|
| Rivara 1984, 1991; Rivara–Levin 1992 | Bisseção pela aresta mais longa (2D/3D/nD) | Nenhuma combinatória; a escolha é geométrica, não "cega" |
| Mitchell 1988/1991 (tese; J. Comput. Appl. Math. 36) | Newest vertex bisection (NVB) em 2D | Aresta de refinamento compatível; em 2D sempre existe (Binev–Dahmen–DeVore 2004 usam emparelhamento perfeito; Biedl et al. 2001 dão o algoritmo) |
| Bänsch 1991; Kossaczký 1994 | NVB em 3D | Kossaczký faz uma pré-subdivisão da malha inicial |
| **Maubach 1995** (SISC 16) | NVB em n dimensões: aresta [v₀, v_γ], tipo cíclico | Malhas "geradas por reflexão" (Kuhn), ou condição de vizinhos refletidos |
| **Traxler 1997** (Computing 59) | Mesmo algoritmo, formulação equivalente | "Reflected domain partition". Para domínio simplesmente conexo: **existe sse toda (n−2)-face interior está num número par de simplexos** |
| Arnold–Mukherjee–Pouly 2000 | Bisseção marcada de tetraedros (5 tipos) | Qualquer malha 3D, mas só refinamentos uniformes de nível múltiplo de n são conformes |
| **Stevenson 2008** (Math. Comp. 77) | Estimativa de fecho (closure) em nD | **Matching neighbour condition** (b), mais fraca que a de Traxler. Não se sabe se toda malha 3D a satisfaz. App. A: pré-subdivisão em (n+1)!/2 simplexos que sempre a satisfaz |
| Alkämper–Gaspoz–Klöfkorn 2018 (SISC 40) | "Compatibilidade fraca": renumeração O(n) de qualquer malha. Garante término, não fecho ótimo | Qualquer |
| Belda-Ferrín, Ruiz-Gironés, Gargallo-Peiró, Roca 2023 (CAD, arXiv 2211.08427) | Bisseção marcada em n estágios, depois "cast" para simplexos refletidos e troca para Maubach | Qualquer malha conforme. Afirmam: "para n > 3 não se conhece método que extraia uma estrutura de reflexão de uma malha arbitrária" |
| Belda-Ferrín et al. IMR 2018 (4D), IMR 2022 (3D, cota ótima de similaridade) | Precursores do anterior | — |
| **Diening–Gehring–Storn 2025** (FoCM, arXiv 2306.02674) | Inicialização de Maubach por **(N+1)-coloração** dos vértices. Algoritmo guloso O(#vértices), N ≤ grau máximo. Provam regularidade de forma e estimativa de fecho de BDV | Qualquer malha. N = n ⇔ condição de Traxler |
| Diening–Storn–Tscherpel 2023 (arXiv 2305.05742) | *Grading* ótimo; geração de vértices definida pela cor, gen(v) = −c(v) | (n+1)-colorida |
| Gehring 2023 (arXiv 2305.03733) | Constante do teorema BDV–Stevenson | — |
| Li 2026 (arXiv 2609.29616) | Primeira estimativa de fecho incondicional para AMP/Bänsch em qualquer malha 3D | Qualquer (3D) |

### 1.2 Como a tese se encaixa (comparação honesta)

**Definição 9 da tese (invariantes de Maubach)** comparada com Stevenson (b):

- (MaubachA) diz que duas células vizinhas do mesmo nível ou têm i = j, ou têm {i, j} = {0, type}. Isso é o análogo, na convenção de Maubach (aresta [0, type]), do "Tᴿ" de Stevenson, que usa a convenção de Traxler (aresta x₀xₙ). Lá, Tᴿ = (xₙ, x₁, …, x_γ, xₙ₋₁, …, x_{γ+1}, x₀)_γ é o simplexo rotulado *com os mesmos filhos* que T, e dois vizinhos são refletidos quando a sequência de T *ou* de Tᴿ coincide com a de T′ em todas as posições menos uma. **[verificar a equivalência exata entre (MaubachA) e "refletidos" de Stevenson, traduzindo as convenções]**
- (MaubachB) cobre vizinhos de níveis adjacentes.
- A diferença estrutural é que **a tese formula um invariante que vale em todos os níveis**, e prova (Lema 7 / Teorema 8) que o algoritmo o preserva. Stevenson e Traxler impõem a condição **só na malha inicial** e depois provam conformidade dos refinamentos uniformes (Teorema 4.3 de Stevenson) e do fecho.

  As duas abordagens provam essencialmente o mesmo teorema. A da tese é local, indutiva e puramente combinatória, no mesmo espírito do "esquema = invariante + algoritmo" que também cobre Mitchell. **Sugestão:** apresentar a Definição 9 como uma forma *hereditária* da matching condition e citar Stevenson (Rem. 4.4: a condição é também necessária) para mostrar que o invariante é o mais fraco possível no nível 0.

**Teoremas 9 e 10 (produto de simplexos; subdivisão baricêntrica):**

- Os dois complexos são **balanceados** no sentido de Stanley: o rank no reticulado e a dimensão da bandeira são (n+1)-colorações. "(n+1)-colorível ⇒ refletido com i = j" é Traxler 1997, reformulado em DGS §3.1.
- Datas: a tese é de 2006, Stevenson 2008 é posterior e Traxler 1997 é anterior. Não dá para reivindicar prioridade sobre o fato geral. Dá para apresentar os Teoremas 9–10 como *instâncias* com provas curtas no formalismo estelar.
- **Ponto de custo favorável à tese?** Não: a pré-subdivisão de Kossaczký/Stevenson usa (n+1)!/2 células contra (n+1)! da baricêntrica, então a deles é 2× menor. Mas DGS mostram que nenhuma das duas é necessária. Para malhas não balanceadas, a comparação justa é com DGS (N > n, sem pré-subdivisão).

**Onde há espaço para contribuição real:**

1. **Variedades fechadas / sementes combinatórias.** Toda a literatura acima trata de domínios em Rⁿ (MEF). Ninguém estudou sementes que são triangulações de variedades fechadas (CP², S²×S²) nem a pergunta "quais triangulações mínimas são balanceadas". Esse é território da combinatória (Klee–Novik, Izmestiev–Klee–Novik, Venturello, Zheng) e ninguém faz a ponte com NVB.
2. **Gaifullin é balanceada**, com 15 vértices. *Correção (2026-09-29):* isso já está no próprio artigo de Gaifullin, que dá a coloração regular dos vértices em 5 cores e prova, pela Proposição 2.3 dele, que X é minimal na classe T_colour(CP²) das triangulações balanceadas. A novidade é só a ponte com a bissecção de Maubach: as sementes fortemente compatíveis são exatamente T_colour, então X é a menor semente possível para CP². Kühnel (9 vértices) não é balanceada.
3. **Esquemas estelares como linguagem unificadora.** Mitchell (faceta) e Maubach (aresta) cabem na mesma forma. Isso conversa com os movimentos estelares (Alexander/Newman; Lickorish 1999), com as gramáticas estelares de Velho (SGP 2003) e 4-8 (Velho–Zorin 2001, que é NVB em 2D) e com os *cross-flips* balanceados (Izmestiev–Klee–Novik 2017), que são o análogo balanceado dos movimentos de Pachner. Frase possível: "a bisseção de Maubach é uma sequência de subdivisões estelares de arestas que preserva balanceamento no sentido fraco da Definição 9".

### 1.3 Critérios de colorabilidade (para a seção sobre malhas iniciais)

- **Joswig 2002** (Math. Z.): um complexo localmente fortemente conexo é balanceado ⇔ seu grupo de projetividades é trivial. Para variedades simplesmente conexas, isso se reduz a toda (n−2)-face ter grau par, o mesmo critério de Traxler.
- Para CP² (simplesmente conexo) basta verificar os graus dos triângulos. Gaifullin: 180 triângulos de grau 4 e 60 de grau 6, todos pares ✔.
- Complexos k-vizinhais com mais de n+1 vértices nunca são balanceados. Isso mata Kühnel imediatamente.

---

## 2. Representações sem ponteiros e busca de vizinhos

| Trabalho | Hierarquia | Semente | Conforme? | Dimensão | Relação com glpt |
|---|---|---|---|---|---|
| Hebert 1994 (J. Symb. Comput. 17) | Bisseção de tetraedros, refinamento "simbólico" | Cubo | sim | 3 | Precursor da codificação simbólica |
| Maubach 1996 (IMR) | Localização de vizinhos em malhas nD de Maubach | Kuhn | sim | n | Precursor direto de Atalay–Mount |
| Evans–Kirkpatrick–Townsend 2001 (Algorithmica 30) | Right-triangulated irregular networks: rótulo binário ⇒ geometria | Quadrado | sim | 2 | Mesma ideia em 2D |
| Lee–De Floriani–Samet 2001 (SMI) | Vizinho em tempo constante, tetraedros de bisseção | Cubo (6 tets) | sim | 3 | Mesma ideia em 3D |
| **Atalay–Mount 2007** (IJCGA 17(6):595–631) | **Código LPT**: geometria e vizinhos por tabela, sem ponteiros | Cubo de Kuhn (órbita de Sym(n)) | sim | n | **Base do glpt** |
| Weiss–De Floriani 2009 (SGP), 2010 (IMR), 2011 (CGF) | Hierarquias de diamantes, uma alternativa de primitiva ao simplexo | Kuhn | sim | n | Visão "diamante = estrela da aresta de bisseção", útil para o fecho |
| Burstedde–Holke 2016 (SISC 38) / t8code | **TM-index**: curva de preenchimento Morton para simplexos (refinamento *red* de Bey); **floresta** sobre malha grossa arbitrária | **Qualquer malha grossa** | **não** (hanging nodes) | 2, 3 | **Análogo arquitetural mais próximo**: id da árvore (≈ `seed_index`) + código linear por árvore + tabela de conectividade entre árvores |
| Burstedde–Wilcox–Ghattas 2011 (p4est) | Floresta de octrees | Qualquer malha de hexaedros | não | 2, 3 | Mesmo padrão "floresta" |
| Bader et al. 2010; Meister–Rahnema–Bader 2016 (sam(oa)²) | NVB 2D com curva de Sierpiński, em pilhas | Malha de triângulos | sim | 2 | Sem ponteiros por ordenação, não por código |
| Dupuy 2020 (PACMCGIT) — Concurrent Binary Trees | Árvore binária implícita num heap de bits; LEB em GPU | Quadrado/malha 2D | sim | 2 | Paralelismo; mostra demanda prática por representações sem ponteiros |
| Boissonnat–Kachanovich–Wintraecken 2021/2023 (SoCG / SICOMP) | Representação permutaédrica de simplexos de Freudenthal–Kuhn e Coxeter; faces e cofaces sensíveis à saída | ℝᵈ inteiro (triangulação uniforme) | sim | d | Codificação compacta de simplexos de Kuhn, mas **sem adaptatividade** |

**Posicionamento do glpt:**
- *Contra Atalay–Mount:* eles dependem de a semente ser uma única órbita de Sym(n) (o cubo). Mostramos que basta uma **semente balanceada**. Cada célula-semente, com os vértices ordenados pela cor, é identificada com o simplexo de referência de Kuhn. O código LPT dentro da semente não muda, e cruzar a fronteira da semente preserva o índice local (i = j), resolvido por uma tabela de adjacência da semente de tamanho 108 × 5 no caso de Gaifullin.

  Isso explica *por que* "a adaptação é fácil". Vale também dizer onde não é: sementes só compatíveis via {0, type} (Stevenson) precisariam de uma correção de permutação na travessia. Sementes DGS com N > n exigiriam códigos de um N-simplexo virtual, e o Remark 12 de DGS avisa que a restrição **não** coincide com refinar o N-simplexo virtual. Isso fica como trabalho futuro.
- *Contra t8code:* a mesma arquitetura (floresta + código linear), mas com bisseção **conforme** (graded, sem hanging nodes) e em dimensão arbitrária (usamos 4).
- *Contra BKW:* eles codificam a triangulação uniforme; o glpt codifica a hierarquia adaptativa.
- Limitação a declarar: o caso `DEEP_CROSSING` (vizinho refinado mais fundo) não está implementado, e a memória é O(folhas) numa tabela hash.

---

## 3. Aplicação: triangulação de superfícies de Riemann / curvas complexas

### 3.1 Computação numérica de superfícies de Riemann (sem malha 2D explícita)
- Deconinck–van Hoeij 2001 (Physica D; pacote Maple `algcurves`) e Deconinck–Patterson (capítulo em Bobenko–Klein, LNM 2013): monodromia, matrizes de período.
- Frauendiener–Klein (Matlab), Bruin–Sijsling–Zotine (Sage `RiemannSurface`), Kranich (tese TUM 2016: continuação homotópica certificada de curvas planas).
- Bobenko–Klein (eds.), *Computational Approach to Riemann Surfaces*, LNM 2013.
- **Comparação:** esses trabalhos tratam a superfície como recobrimento ramificado de ℂP¹ e calculam invariantes analíticos. Nós construímos um **complexo simplicial** da curva em CP², o que é complementar (daria para usar nossa malha para verificar gênero via χ).

### 3.2 Visualização de superfícies de Riemann
- Hanson 1994 (Notices AMS 41): curvas de Fermat em CP² projetadas em 3D, **paramétricas** e explorando simetria.
- Trott (Mathematica), Nieser–Poelke–Polthier 2010 (GMP, LNCS 6130): malha como gráfico multivalorado sobre ℂ, a partir dos pontos de ramificação, com cortes.
- Kranich 2015 (arXiv 1507.04571): superfícies de Riemann coloridas por domínio em GPU.
- **Comparação:** todos usam uma projeção para ℂ (cartas, cortes). Nós somos **intrínsecos em CP²**: o infinito não é especial e não há cortes. O preço é o problema de visualização por cartas (ver memória *viz chart consistency*).

### 3.3 Variedades implícitas de codimensão > 1 (a família técnica mais próxima)
- **Allgower–Schmidt 1985** (SINUM 22): aproximação PL de variedades implícitas de codimensão arbitrária sobre Freudenthal–Kuhn, com refinamento. Allgower–Georg (livro, SIAM Classics 2003): PL-continuação e pivoteamento. Todd J₃, Eaves–Saigal: triangulações com refinamento para homotopias simpliciais. **Essa é a linhagem histórica do nosso "continuação + Kuhn".**
- **Weigle–Banks 1996** (IEEE Vis): "complex-valued contour meshing": o zero de f: ℂ² → ℂ como superfície em ℝ⁴, por extração simplicial em grade. **O antecessor mais direto da nossa extração de codimensão 2**, mas em grade afim uniforme de ℝ⁴, não em CP².
- Bhaniramka–Wenger–Crawfis 2004 (TVCG 10): isosuperfícies em qualquer dimensão (codimensão 1), com prova de variedade.
- Castelo, Nonato, Siqueira, Minghim, Tavares 2006 (C&G 30): triangulação J1a, adaptativa em qualquer dimensão, para variedades implícitas. Linhagem brasileira (Tavares, também em Gomes–Tavares 1989), concorrente natural na adaptatividade nD.
- **Boissonnat–Kachanovich–Wintraecken 2021/2023:** traçado de isovariedades a partir de um ponto semente (continuação!) em Freudenthal–Kuhn/Coxeter, em tempo polinomial em d, com garantias. Boissonnat–Wintraecken 2022 (FoCM): corretude topológica de aproximações PL de isovariedades. **Referência obrigatória:** mesma estratégia (seguir a variedade de célula vizinha em célula vizinha), triangulação uniforme, com garantias. Nós: adaptativo, em CP², sem garantias.
- Ju, Du, Zhou, Carr, Ju 2024 (SIGGRAPH/TOG): grade simplicial adaptativa para complexos implícitos, incluindo **redes de curvas** (codimensão 2 em ℝ³), a partir de uma semente cubo-de-6-tets. **[verificar a regra exata de bisseção e o critério]**
- Gerstner–Rumpf 2000; Gregorski et al. 2002 (Vis); Pascucci 2002 (SGS, ⁿ√2): extração de isosuperfícies sobre hierarquias de bisseção. Mostram que a combinação "hierarquia de Maubach + extração" já é usada em visualização, mas só em codimensão 1 e ℝ³.

### 3.4 Certificação (para situar a taxa de células ruins)
- Plantinga–Vegter 2004 (SGP): isotopia garantida via aritmética intervalar em quadtree/octree.
- Mourrain–Pavone 2009 (JSC 44): subdivisão em base de Bernstein. Reuter et al. 2008 (Vis. Comput. 24): **Bernstein baricêntrico em simplexos**, que é o paralelo exato do nosso `--bernstein` sobre 2-faces.
- Burr–Choi–Galehouse–Yap 2012 (DCG); Alberti–Mourrain–Wintz 2008 (CAGD): topologia de curvas implícitas por subdivisão.
- Bertini_real (Brake et al. 2017, TOMS 44): decomposição celular numérica de curvas e superfícies reais. Uma curva complexa é uma superfície real em ℝ⁴, então é um concorrente possível **[verificar se lida com codimensão 2 em ℝ⁴]**.
- **Mensagem para o artigo:** nosso método é heurístico e adaptativo. `--tangency-threshold` e `--repair-rounds` reduzem a taxa de células ruins (melhor resultado: 0,180 %), mas não a zeram. Posicionar como "triangulação prática, verificável a posteriori (χ = 2 − 2g)", não como certificada.

---

## 4. Triangulações de CP² (para a seção sobre a semente)
- Kühnel–Banchoff 1983 (Math. Intelligencer 5): CP²₉, único (Kühnel–Lassmann; prova sem computador de Arnoux–Marin). 3-vizinhal ⇒ não balanceada.
- Banchoff–Kühnel 1992 (Geom. Dedicata): "equilibrium triangulations" de CP². **[verificar relevância]**
- **Gaifullin 2009** (Trudy MIAN 266 / Proc. Steklov Inst. 266:29–48; arXiv 0904.4222): *"A minimal triangulation of complex projective plane admitting a chess colouring of four-dimensional simplices"*. 15 vértices, automorfismos S₄ × S₃ agindo por **isometrias de Fubini–Study**. Isso explica a razão de arestas FS = 1,4636 com variância zero que medimos. Traz também uma subdivisão com 33 vértices em que a aplicação momento é simplicial. ⚠ **O relatório atual cita arXiv:0904.4222 com título errado** ("Construction of combinatorial manifolds with prescribed sets of links of vertices" é outro artigo de Gaifullin). Corrigir.
- Schwartz 2022 (Math. Intelligencer): trissecção do CP²₉.
- Effenberger–Spreer (simpcomp), GAP: ferramentas de verificação que já usamos.

---

## 5. Sugestões concretas para o artigo

1. **Enunciado central:** "Maubach bisection works verbatim on any *balanced* (n+1)-colorable initial triangulation, including triangulations of closed manifolds. Pointerless LPT codes extend to such seeds by adding a seed index and a seed adjacency table." Deixar claro que a primeira metade é conhecida (Traxler, DGS) e que as contribuições são:
   - o formalismo de esquemas estelares com invariante hereditário;
   - sementes em variedades fechadas (CP² de Gaifullin, balanceada);
   - a extensão do LPT;
   - a aplicação.
2. **Experimentos de comparação que valem a pena:**
   - Kühnel com DGS (N = 8) contra Gaifullin (N = 4): células, qualidade FS, taxa de células ruins. Seria o primeiro uso de DGS fora de ℝⁿ.
   - Barycêntrica(Kühnel) contra a pré-subdivisão de Stevenson (Kühnel) contra Gaifullin: custo em células e qualidade.
   - Weigle–Banks (grade de ℝ⁴ numa carta afim) contra CP² intrínseco, na mesma curva: células visitadas e artefatos perto do infinito.
3. **Correções no relatório existente** antes de reaproveitar texto:
   - Gaifullin *tem* estrutura de rank global (uma coloração);
   - título e referência de Gaifullin;
   - a §3 "coloring problem on the dual graph" deve ser reinterpretada como busca de (n+1)-coloração dos vértices.
4. **Pergunta aberta atraente:** menor triangulação balanceada de CP²; e se a subdivisão com 33 vértices de Gaifullin (aplicação momento simplicial) dá uma semente ainda melhor.
