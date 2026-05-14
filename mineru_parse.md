## 选择区域 1

**来源页面:** [输出文件格式 - MinerU](https://opendatalab.github.io/MinerU/zh/reference/output_files/)

**选择器信息:**
- XPath: `/html/body/div[3]/main[1]/div[1]/div[3]/article[1]/div[1]`
- CSS Selector: `html.js-focus-visible.js > body > div.md-container:nth-of-type(3) > main.md-main > div.md-main__inner.md-grid > div.md-content:nth-of-type(3) > article.md-content__inner.md-typeset > div.md-selector-preview.md-selector-highlight`

# MinerU 输出文件说明
 ## 概览
 `mineru` 命令执行后，除了输出主要的 markdown 文件外，还会生成多个辅助文件用于调试、质检和进一步处理。这些文件包括：
 具体会生成哪些文件，取决于后端类型和输入文档类型。
 - **可视化调试文件**：帮助用户直观了解文档解析过程和结果
- **结构化数据文件**：包含详细的解析数据，可用于二次开发

 ### 通用内容列表 V2 (content_list_v2.json)(开发中，格式可能调整)
 **文件命名格式**： `{原文件名}_content_list_v2.json`
 ##### 功能说明
 `content_list_v2.json` 是 3.0 起新增的结构化输出文件，所有后端都会在保留 `content_list.json` 的同时额外输出该文件：
 - 顶层是按页分组的列表，便于按页消费结果
- 每个内容块使用统一的 `type + content` 结构，适合程序化处理
- 不同后端和输入类型支持的 `type` 会有所不同
 ##### 通用字段
 字段名 类型 说明 `type` `string` 内容类型 `content` `dict` 与 `type` 对应的结构化内容 `bbox` `list[int]` 可选，0-1000 范围的边界框 `anchor` `string` 可选，部分 `DOCX` 标题或索引项会携带锚点
 其中 `image` / `chart` 类型还可能包含可选顶层字段 `sub_type`，用于表示视觉子类型。
 ##### 常见类型
 类型 说明 `title` 标题块，包含 `title_content` 与 `level` `paragraph` 段落块，包含 `paragraph_content` `equation_interline` 行间公式，包含 `math_content`、 `math_type` `image` / `table` / `chart` / `seal` 视觉类块，包含图片路径、说明文字等结构化字段 `code` 代码块，包含 `code_content`、 `code_caption`、 `code_footnote`、 `code_language` `algorithm` 算法块，包含 `algorithm_content`、 `algorithm_caption`、 `algorithm_footnote` `list` / `index` 列表与索引，包含 `list_items` `page_header` / `page_footer` / `page_number` / `page_aside_text` / `page_footnote` 页面辅助块
 ##### 示例数据
 ```
[
    [
        {
            "type": "title",
            "content": {
                "title_content": [
                    {
                        "type": "text",
                        "content": "1 Introduction"
                    }
                ],
                "level": 1
            },
            "bbox": [
                83,
                121,
                917,
                156
            ]
        },
        {
            "type": "page_footnote",
            "content": {
                "page_footnote_content": [
                    {
                        "type": "text",
                        "content": "* Corresponding author"
                    }
                ]
            },
            "bbox": [
                71,
                815,
                915,
                841
            ]
        }
    ]
]

```
 