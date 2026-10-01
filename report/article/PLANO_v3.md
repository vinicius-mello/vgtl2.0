# Plano: artigo v3, centrado na triangulação de superfícies de Riemann em CP²

*2026-10-01.* A versão anterior está guardada na tag `article-v2` (commit c709a4b) e em `article_v2.pdf`. A v1 está na tag `article-v1` e em `article_v1.pdf`.

## 1. Mudança de eixo

| | v2 (atual) | v3 (proposta) |
|---|---|---|
| Pergunta | Por que Maubach está correto e em quais sementes ele roda? | Como triangular uma curva plana suave **intrinsecamente** em CP²? |
| Contribuição | Teoria (esquemas estelares, sementes balanceadas) + LPT | Um **método** completo e medido, montado com técnicas conhecidas |
| Aplicação | Exemplo curto (§7) com medições no apêndice | É o artigo |
| Teoria de Maubach | Provas completas (§3–4) | Citada (Maubach, Traxler, Stevenson, DGS), sem provas |

Mensagem central proposta: *uma curva plana suave é uma superfície real compacta em CP² e pode ser triangulada sem cartas, sem cortes de ramificação e sem tratamento especial do infinito. Basta refinar adaptativamente uma triangulação de CP² e seguir a curva de célula vizinha em célula vizinha. Todas as peças existem; a combinação é nova, e o resultado pode ser verificado a posteriori (χ = 2 − 2g).*

Título candidato: *Intrinsic triangulation of plane algebraic curves in CP² by adaptive bisection*, ou *Triangulating Riemann surfaces inside CP²*.

## 2. As peças e de onde vêm

| Peça | Origem (citar) | O que é nosso |
|---|---|---|
| Semente de CP² com 15 vértices | Gaifullin 2009 (coordenadas, coloração, minimalidade) | escolha e rotação unitária aleatória |
| Compatibilidade da semente | Traxler 1997 / DGS 2025: uma semente balanceada é compatível | só a observação (1 parágrafo) |
| Bisseção conforme | Maubach 1995, Traxler 1997, Stevenson 2008 | pontos médios em geodésicas de Fubini–Study |
| Malha sem ponteiros | Atalay–Mount 2007 (LPT) | extensão a semente não cúbica (ver decisão D1) |
| Extração de codimensão 2 | Weigle–Banks 1996; Allgower–Schmidt 1985 | Newton por 2-face na carta mais bem condicionada, pareamento por comprimento geodésico |
| Continuação | Allgower–Georg; Boissonnat–Kachanovich–Wintraecken | continuação por vizinhança combinatória na malha adaptativa |
| Refino perto de tangências | — (heurísticas) | escore de transversalidade (Cauchy–Riemann), rodadas de reparo |
| Verificação topológica | fórmula de Euler | divisão de vértices pinçados, fechamento por cones com Newton |
| Visualização | Hanson 1994; Dutter 2026 (mapa sem cartas) | cartas planas, `render_obj.py` |

## 3. O que sai, o que fica, o que muda

**Sai (fica preservado na tag `article-v2`):**
- §3 inteira (esquemas de subdivisão estelar, Mitchell) e as provas da §4 (lemas de direção, filhos e duas células, Teorema de Maubach, algoritmo, relação com Stevenson).
- §5.1–5.3 como seção própria: a equivalência compatível ⇔ balanceada e as construções clássicas viram **um parágrafo** com citações.
- Contribuições 1 e 2 da introdução atual.

**Fica, condensado:**
- §5.4 *Seeds for CP²*, reduzida a uma subseção "A semente":
  - Gaifullin é balanceada e, portanto, compatível;
  - Kühnel não é balanceada (é 3-vizinhal);
  - medidas de qualidade (razão de arestas FS 1,4636 com variância zero; 108 células).
- Observação "index-only matching não basta": é um alerta prático útil (foi o nosso erro com Kühnel). Pode virar nota de rodapé ou ir para o apêndice.

**Fica e cresce:** §7 (Exemplo) e o Apêndice A viram o corpo do artigo.

## 4. Estrutura proposta para a v3

1. **Introdução**: o problema, por que fazer intrinsecamente (infinito, cortes, cúspide de canto na compactificação S²×S², ver item E3), contribuições, figura principal.
2. **Trabalhos relacionados**: reaproveitar a §3 de `related_work.md` (superfícies de Riemann numéricas, visualização, variedades implícitas de codimensão 2, certificação). A parte de bisseção encolhe para um parágrafo.
3. **Preliminares**: CP², métrica de Fubini–Study, geodésicas e pontos médios; curvas planas e gênero; bisseção de Maubach (enunciado e citações, sem prova).
4. **Malha adaptativa de CP²**: semente de Gaifullin, compatibilidade via coloração, rotação aleatória, bisseção ao longo de geodésicas FS, malha sem ponteiros (conforme D1).
5. **Extração e continuação**: cruzamentos nas 2-faces, arcos nas facetas, polígonos, células ruins, continuação, limiar de tangência, rodadas de reparo.
6. **Verificação topológica e fechamento**: vértices pinçados, estimativa do gênero, fechamento por cones.
7. **Experimentos**:
   - taxa de células ruins × profundidade;
   - causas das falhas;
   - robustez a rotações;
   - gênero bruto e após fechamento;
   - precisão dos vértices;
   - tempo e memória;
   - itens novos da seção 5.
8. **Visualização**: cartas planas e o mapa de Dutter (figuras atuais).
9. **Limitações e conclusão**: sem certificação, curvas singulares, componentes pequenas perdidas, grau alto.

## 5. Experimentos

**Já existem** (dados e scripts em `figures/`):
- `badrate.csv` / `badrate.pdf`;
- `rotations.csv`;
- tabelas de gênero, causas e precisão do Apêndice A;
- figuras `elliptic_flat` e `quartic_alpha`.

**Novos, por ordem de valor para um artigo de método:**
- **E1. Comparação com a abordagem por carta afim (estilo Weigle–Banks).** Mesma curva, grade simplicial de ℝ⁴ numa carta, contra CP² intrínseco. Medir células visitadas, artefatos perto do infinito e o que se perde fora da janela. É a comparação que um revisor vai pedir.
- **E2. Tempo e memória.** Hoje há memória (nmt 222 contra glpt 19,5 B/folha), mas falta uma tabela de tempo por profundidade e por curva.
- **E3. S²×S² contra CP².** `examples/top/riemann/riemann.cpp` já triangula em S²×S². Lá a curva elíptica tem uma cúspide genuína em (∞,∞), e em CP² não. É um bom argumento motivador, com dados que já sabemos produzir.
- **E4. Mais curvas.** Grau 5–6, coeficientes aleatórios, gênero até 10. Mostra quando a taxa de células ruins começa a pesar.
- **E5 (opcional). Kühnel com inicialização DGS (N = 8) contra Gaifullin.** É um teste de semente, interessante mas periférico na nova linha. Talvez caiba numa observação.

## 5b. Desenho do E1 (experimento principal)

**Linha de base.** É a abordagem "Weigle–Banks": triangular {F(x,y,1)=0} numa carta afim ℂ² = ℝ⁴. Para a comparação ser justa, usa-se **o mesmo código** e só muda a geometria:
- **Semente:** caixa [−R,R]⁴ em (Re x, Im x, Re y, Im y), triangulada por Kuhn. O bloco refletido 2⁴ tem 384 células e é balanceado, então o glpt roda sem alteração; falta só a tabela de adjacência da semente.
- **Pontos:** os vértices são (x, y, 1) ∈ ℂ³, e os pontos médios são afins na carta, não de Fubini–Study.
- **Resto do pipeline idêntico:** extração por 2-face, continuação, tangência, reparo e `--bernstein`.

Implementação: opções `--seed box:R` e `--midpoint affine` no `riemann_cp2_glpt`.

**Métricas.**
1. **Cobertura**, pela área de Fubini–Study da superfície extraída. Pelo teorema de Wirtinger, uma curva de grau n tem área FS exatamente n·área(reta), então a cobertura é área / (n·área(reta)). A constante se calibra rodando uma reta (grau 1). No método intrínseco isso também **valida** a extração, porque deve dar ≈ 1. Na caixa, mede quanto da curva se perde fora dela, em função de R.
2. **Topologia:**
   - em ℂP², a superfície é fechada e χ = 2 − 2g;
   - na caixa, a superfície tem bordo e o gênero não se recupera diretamente. Contar os laços de bordo e compará-los com os n pontos no infinito (F(x,y,0)=0).
3. **Custo para a mesma resolução.** Refinar até que o diâmetro FS das células que tocam a curva fique < h, e não até uma profundidade fixa. A grade afim desperdiça células perto da origem e é grossa, em FS, perto da borda da caixa. Comparar folhas, tempo e taxa de células ruins para os mesmos valores de h.
4. **Distorção:** histograma da razão de arestas FS das células extraídas, nos dois métodos.

**Variações.** R ∈ {2, 4, 8, 16} (cobertura → 1 devagar e custo crescendo), curvas elíptica e quártica de Fermat, e 1–2 rotações. Mencionar, sem implementar, que cobrir ℂP² com três cartas exige costurar as sobreposições, que é exatamente o que o método intrínseco evita.

## 6. Decisões

- **D1 (decidido 2026-10-01). A extensão do LPT a sementes balanceadas fica**, como contribuição própria, na seção 4 (malha adaptativa de CP²). Inclui:
  - o teste de fronteira por tabela;
  - a travessia de semente com o mesmo código;
  - a busca de vizinhos totalmente combinatória;
  - a verificação exata (5 milhões de consultas);
  - a comparação de memória (nmt 222 contra glpt 19,5 B/folha).

  É a única peça nova do lado da malha, e a introdução deve listá-la entre as contribuições.
- **D2 (decidido 2026-10-01). Ênfase em E1**: a comparação com a abordagem por carta afim é o experimento principal, e E4 fica secundário. O público é de geometria computacional / computação gráfica (C&G, CAGD, SMI, CGF); o veículo exato ainda não foi escolhido.
- **D3 (revisto 2026-10-01). `--bernstein` sai do artigo**: só funciona com faces planas na carta, que deixaram de ser o padrão; a certificação por aritmética intervalar nas faces geodésicas vira trabalho futuro. (Decisão anterior: `--bernstein` fica, com o experimento E6.)
- **D4 (decidido 2026-10-01). A teoria é arquivada para outro momento**, na tag `article-v2` e em `article_v2.pdf`. Na v3 ela não aparece nem como trabalho futuro. A nota `report/face_scheme/` fica como registro.

**E6 (cancelado com a revisão do D3).** No `riemann_cp2_glpt`, comparar `--bernstein` com o modo padrão:
- taxa de células ruins, folhas, tempo e memória;
- profundidades 14–20, com e sem tangência e reparo;
- curvas elíptica e quártica.

As medidas antigas com `--bernstein` (5,78 %, 3,20 % e 2,24 % nas profundidades 10, 14 e 16) são do `riemann_cp2` com nmt e não valem diretamente para o glpt com continuação.

## 7. Ordem de trabalho sugerida

1. Decisões D1–D4 tomadas; falta só escolher o veículo exato.
2. Esqueleto da v3 em `article.tex`: nova estrutura, mover o texto reaproveitável da v2 (§5.4, §6 condensada, §7, Apêndice A) e cortar §3–4.
3. Reescrever a introdução e os trabalhos relacionados para o novo eixo.
4. Implementar e rodar E1 (prioridade; exige a semente em caixa e os pontos médios afins), depois E2, E3 e E6 (baratos, o código já existe) e por fim E4.
5. Revisar figuras e tabelas, conclusão e créditos.
6. Atualizar `related_work.md` (seção 5, "Sugestões") para refletir a nova linha.
